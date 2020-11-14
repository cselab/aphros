// Created by Petr Karnakov on 11.11.2020
// Copyright 2020 ETH Zurich

#include <algorithm>
#include <iostream>
#include <map>

#include "distr/distrbasic.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/format.h"
#include "util/timer.h"

#include "commmap.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

std::string Gather(const std::string& str, MPI_Comm comm) {
  const std::vector<char> buf(str.begin(), str.end());
  const int commsize = MpiWrapper::GetCommSize(comm);
  int bufsize = buf.size();
  if (MpiWrapper::IsRoot(comm)) {
    std::vector<int> sizes(commsize);

    MPI_Gather(&bufsize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> offsets = {0};
    for (auto& q : sizes) {
      offsets.push_back(offsets.back() + q);
    }
    std::vector<char> bufall(offsets.back());

    MPI_Gatherv(
        buf.data(), buf.size(), MPI_CHAR, bufall.data(), sizes.data(),
        offsets.data(), MPI_CHAR, 0, comm);
    return {bufall.begin(), bufall.end()};
  } else {
    MPI_Gather(&bufsize, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);

    MPI_Gatherv(
        buf.data(), buf.size(), MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR,
        0, comm);
    return "";
  }
}

class StreamMpi {
 public:
  StreamMpi(std::ostream& out, MPI_Comm comm) : out_(out), comm_(comm) {}
  ~StreamMpi() {
    Flush();
  }
  void Flush() {
    auto str = Gather(buf_.str(), comm_);
    if (MpiWrapper::IsRoot(comm_)) {
      out_ << str;
    }
    buf_ = {};
  }
  template <class T>
  StreamMpi& operator<<(const T& value) {
    buf_ << value;
    return *this;
  }

 private:
  std::ostream& out_;
  MPI_Comm comm_;
  std::stringstream buf_;
};

std::vector<std::string> FieldToStrings(
    const FieldCell<Scal>& fcu, const M& m, std::string fmt = "{:4.1f}") {
  auto wb = m.GetSuBlockCells().GetBegin();
  auto we = m.GetSuBlockCells().GetEnd();
  std::vector<std::string> res;
  for (int y = we[1]; y-- > wb[1];) {
    if (y == wb[1] || y + 2 == we[1]) {
      res.emplace_back("----");
      const int sepwidth = util::Format(fmt, "0").length();
      for (int x = 0; x < (we[0] - wb[0]) * sepwidth; ++x) {
        res.back() += "-";
      }
    }
    res.emplace_back();
    for (int x = wb[0]; x < we[0]; ++x) {
      const IdxCell c = m.GetIndexCells().GetIdx(MIdx(x, y, 0));
      const bool xedge = x == wb[0] || x + 2 == we[0];
      res.back() += util::Format(fmt + "{}", fcu[c], xedge ? " |" : "");
    }
  }
  return res;
}

void PrintField(
    std::ostream& out, const FieldCell<Scal>& fcu, M& m,
    std::string fmt = "{:4.1f}") {
  auto sem = m.GetSem(__func__);
  struct {
    std::vector<std::string> res;
    std::vector<std::vector<char>> resbuf;
    std::vector<int> meshid;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    t.res = FieldToStrings(fcu, m, fmt);
    for (auto s : t.res) {
      t.resbuf.emplace_back(s.begin(), s.end());
    }
    t.meshid.push_back(m.GetId());
    m.Reduce(&t.resbuf, Reduction::concat);
    m.Reduce(&t.meshid, Reduction::concat);
  }
  if (sem()) {
    if (m.IsRoot()) {
      const int nlines = t.res.size();
      std::map<int, size_t> meshid_to_offset;
      for (size_t i = 0; i < t.meshid.size(); ++i) {
        meshid_to_offset[t.meshid[i]] = i;
      }
      const MIdx blocks = m.flags.global_blocks;
      for (int y = blocks[1]; y-- > 0;) {
        for (int i = 0; i < nlines; ++i) {
          bool first = true;
          for (int x = 0; x < blocks[0]; ++x) {
            auto offset = meshid_to_offset[blocks[0] * y + x];
            auto& sbuf = t.resbuf[nlines * offset + i];
            if (!first) {
              out << "  ";
            } else {
              first = false;
            }
            out << std::string(sbuf.begin(), sbuf.end());
          }
          out << "\n";
        }
        out << "\n";
      }
    }
  }
}

void Run(M& m, Vars&) {
  auto sem = m.GetSem(__func__);
  struct {
    FieldCell<Scal> fc_rank;
    FieldCell<Scal> fc_flat;
    FieldCell<Scal> fc_nci;
    FieldCell<bool> fc_has_neighbors;
    std::vector<IdxCell> flat_to_cell;
    std::vector<std::array<char, 6>> flat_to_nci_inner;
    using Expr = typename M::Expr;
    FieldCell<Expr> fc_system;
  } * ctx(sem);

  static int flat_current = 0;
  static std::unique_ptr<StreamMpi> stream_mpi;
  static std::vector<Scal> amgx_data;
  static std::vector<int> amgx_cols;
  static std::vector<int> amgx_row_ptrs;
  static std::map<int, std::vector<int>> amgx_recv;
  static std::map<int, std::vector<int>> amgx_send;

  auto& t = *ctx;

  auto next_inner = [&t, &m](IdxCell c) {
    t.fc_flat[c] = flat_current++;
    t.flat_to_cell.push_back(c);
    t.flat_to_nci_inner.emplace_back();
    auto& nci = t.flat_to_nci_inner.back();
    size_t j = 0;
    for (auto q : m.Nci(c)) {
      const auto cn = m.GetCell(c, q);
      if (t.fc_has_neighbors[cn] == 0) {
        nci[j++] = q;
      }
    }
    for (auto q : m.Nci(c)) {
      const auto cn = m.GetCell(c, q);
      if (t.fc_has_neighbors[cn]) {
        nci[j++] = q;
      }
    }
    fassert_equal(j, nci.size());
  };

  const auto comm = m.GetMpiComm();
  const int rank = MpiWrapper::GetCommRank(comm);
  if (sem()) {
    if (m.IsLead()) {
      stream_mpi = std::make_unique<StreamMpi>(std::cout, comm);
    }
    t.fc_rank.Reinit(m, 0);
    for (auto c : m.Cells()) {
      t.fc_rank[c] = rank;
    }
    m.Comm(&t.fc_rank);
  }
  if (sem()) {
    t.fc_has_neighbors.Reinit(m, false);
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        const auto cn = m.GetCell(c, q);
        if (t.fc_rank[cn] != rank) {
          t.fc_has_neighbors[c] = true;
        }
      }
    }
  }
  if (sem("flat-inner-no-neghbors")) {
    if (m.IsLead()) {
      // XXX must be initialized so the function is re-entrant
      flat_current = 0;
    }
    t.fc_flat.Reinit(m, -1);
    for (auto c : m.Cells()) {
      // cells that do not have neighbors from remote ranks
      if (!t.fc_has_neighbors[c]) {
        next_inner(c);
      }
    }
  }
  if (sem("flat-inner-other")) {
    for (auto c : m.Cells()) {
      // cells that have a neighbor from remote rank
      if (t.fc_has_neighbors[c]) {
        next_inner(c);
      }
    }
    // fill halo cells with already computed indices
    m.Comm(&t.fc_flat);
  }
  if (sem("flat-halo")) {
    for (auto c : m.Cells()) {
      // halo cells from remote ranks
      // -1 to prevent traversing the same cell twice
      for (auto q : m.Nci(c)) {
        const auto cn = m.GetCell(c, q);
        if (t.fc_rank[cn] != rank) {
          t.fc_flat[cn] = -1;
        }
      }
    }
    for (auto c : m.Cells()) {
      // halo cells from remote ranks
      for (auto q : m.Nci(c)) {
        const auto cn = m.GetCell(c, q);
        if (t.fc_rank[cn] != rank) {
          if (t.fc_flat[cn] == -1) {
            t.fc_flat[cn] = flat_current++;
            t.flat_to_cell.push_back(cn);
          }
        }
      }
    }
  }
  if (sem()) {
    t.fc_nci.Reinit(m, 0);
    for (size_t i = 0; i < t.flat_to_cell.size(); ++i) {
      const auto c = t.flat_to_cell[i];
      Scal w = 0;
      for (auto q : t.flat_to_nci_inner[i]) {
        w = w * 10 + q;
      }
      t.fc_nci[c] = w;
    }

    if (m.IsLead()) {
      amgx_data.clear();
      amgx_cols.clear();
      amgx_row_ptrs = {0};
      amgx_recv.clear();
      amgx_send.clear();
    }

    /*
    t.fc_system.Reinit(m);
    for (size_t i = 0; i < t.flat_to_nci_inner.size(); ++i) {
      const auto c = t.flat_to_cell[i];
      amgx_cols.push_back(t.fc_flat[c]);
      amgx_data.push_back(t.fc_system[c][0]);
      for (auto q : t.flat_to_nci_inner[i]) {
        const auto cn = m.GetCell(c, q);
        amgx_cols.push_back(t.fc_flat[cn]);
        amgx_data.push_back(t.fc_system[c][1 + q]);
      }
      amgx_row_ptrs.push_back(amgx_cols.size());
    }
    */

    // collect halo cells to receive and inner cells to send
    for (size_t i = t.flat_to_nci_inner.size(); i < t.flat_to_cell.size();
         ++i) {
      const auto c = t.flat_to_cell[i];
      fassert(t.fc_rank[c] != rank);
      amgx_recv[t.fc_rank[c]].push_back(i);
      int cnt = 0;
      for (auto q : m.Nci(c)) {
        const auto cn = m.GetCell(c, q);
        if (t.fc_rank[cn] == rank) {
          amgx_send[t.fc_rank[c]].push_back(i);
          ++cnt;
        }
      }
      fassert_equal(cnt, 1);
    }
  }
  if (sem()) {
    if (m.IsLead()) {
      std::vector<int> amgx_neighbors;
      std::vector<int> amgx_recv_sizes;
      std::vector<const int*> amgx_recv_maps;
      for (const auto& p : amgx_recv) {
        amgx_neighbors.push_back(p.first);
        amgx_recv_sizes.push_back(p.second.size());
        amgx_recv_maps.push_back(p.second.data());
      }
      auto& out = *stream_mpi;
      out << "\nrank=" << rank << '\n';
      out << "size=" << t.flat_to_nci_inner.size() << '\n';
      for (const auto& p : amgx_recv) {
        out << "recv from rank " << p.first << '\n';
        for (auto i : p.second) {
          out << i << " ";
        }
        out << '\n';
      }
      for (const auto& p : amgx_send) {
        out << "send to rank " << p.first << '\n';
        for (auto i : p.second) {
          out << i << " ";
        }
        out << '\n';
      }
    }
  }
  if (sem.Nested()) {
    PrintField(std::cout, t.fc_rank, m, "{:2g}");
  }
  if (sem.Nested()) {
    PrintField(std::cout, t.fc_flat, m, "{:3g}");
  }
  if (sem.Nested()) {
    // PrintField(std::cout, t.fc_nci, m, " {:06d}");
  }
  if (sem.Nested()) {
    // PrintField(std::cout, t.fc_has_neighbors, m, "{:2g}");
  }
  if (sem()) {
    if (m.IsLead()) {
      stream_mpi.reset();
    }
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Test for communication maps.", mpi.IsRoot());
  parser.AddVariable<int>("--nx", 16).Help("Mesh size in x-direction");
  parser.AddVariable<int>("--ny", 16).Help("Mesh size in y-direction");
  parser.AddVariable<int>("--block", 8)
      .Help("Block size in x- and y-directions")
      .Options({8, 16, 32});
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const int nx = args.Int["nx"];
  const int ny = args.Int["ny"];
  const int block = args.Int["block"];
  Subdomains<MIdx> sub(
      MIdx(nx, ny, 1), MIdx(block, block, 1), mpi.GetCommSize());
  std::string conf = "";
  conf += sub.GetConfig();
  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
