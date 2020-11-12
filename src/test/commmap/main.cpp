// Created by Petr Karnakov on 11.11.2020
// Copyright 2020 ETH Zurich

#include <iostream>

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
          for (int x = 0; x < blocks[0]; ++x) {
            auto offset = meshid_to_offset[blocks[0] * y + x];
            auto& sbuf = t.resbuf[nlines * offset + i];
            out << std::string(sbuf.begin(), sbuf.end()) << "  ";
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
    FieldCell<Scal> fc_idx;
    FieldCell<Scal> fc_num_neighbors;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    t.fc_rank.Reinit(m, 0);
    for (auto c : m.SuCells()) {
      t.fc_rank[c] = MpiWrapper::GetCommRank(m.GetMpiComm());
    }
    m.Comm(&t.fc_rank);
  }
  if (sem()) {
    t.fc_num_neighbors.Reinit(m, 0);
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        const auto cn = m.GetCell(c, q);
        if (t.fc_rank[c] != t.fc_rank[cn]) {
          ++t.fc_num_neighbors[c];
        }
      }
    }
  }
  static int colidx = 0;
  if (sem()) {
    t.fc_idx.Reinit(m, 0);
    for (auto c : m.Cells()) {
      if (t.fc_num_neighbors[c] == 0) {
        t.fc_idx[c] = colidx++;
      }
    }
  }
  if (sem()) {
    for (auto c : m.Cells()) {
      if (t.fc_num_neighbors[c]) {
        t.fc_idx[c] = colidx++;
      }
    }
  }
  if (sem.Nested()) {
    PrintField(std::cout, t.fc_rank, m, "{:2g}");
  }
  if (sem.Nested()) {
    PrintField(std::cout, t.fc_idx, m, "{:3g}");
  }
  if (sem.Nested()) {
    //PrintField(std::cout, t.fc_num_neighbors, m, "{:2g}");
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
