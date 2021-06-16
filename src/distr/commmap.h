// Created by Petr Karnakov on 15.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include "geom/mesh.h"
#include "util/distr.h"
#include "util/format.h"

/*

Things AMGX needs:

* indices of inner cells that do not refer to neighboring ranks
* indices of inner cells that refer to neighboring ranks
* indices of halo cells and the corresponding ranks to receive
* indices of inner cells and the corresponding ranks to send
* coefficients from inner cells that do not refer to neighboring ranks
* coefficients from inner cells that refer to neighboring ranks

Having the ordering of neighbors Nci() depend on the cell (to position the halo
cells after the inner cells as required by AMGX) is not feasible.
So the reordering of coefficients must be done as a post-processing step.
Can be done separately for each equation, say the array is appended by 7
elements at once and filled in correct ordering.

Another option is to store a sequence of elements from `Nci()` for each cell
(`char` to save memory) in correct ordering.

One way to obtain the correct ordering is to form the system in any order
and then sort it with a stable sort by:
number of neighbors, cell index, neighbor rank.
Needs to be done only once to obtain the indexing.

Mappings needed:
* from IdxCell to row index
* from <IdxCell,nci> to column index
* row_to_cell from row index to IdxCell


*/

template <class It>
std::vector<size_t> ArgSort(It begin, It end) {
  std::vector<size_t> res(end - begin);
  std::iota(res.begin(), res.end(), 0);
  std::stable_sort(res.begin(), res.end(), [begin](auto a, auto b) {
    return *(begin + a) < *(begin + b);
  });
  return res;
}

template <class It>
void Reorder(It begin, const std::vector<size_t>& indices) {
  std::vector<typename decltype(begin)::value_type> values(
      indices.size(), *begin);
  for (size_t i = 0; i < indices.size(); ++i) {
    values[i] = *(begin + indices[i]);
  }
  for (size_t i = 0; i < indices.size(); ++i) {
    *(begin + i) = values[i];
  }
}

template <class It1, class It2>
void SortByFirst(It1 begin1, It1 end1, It2 begin2) {
  auto indices = ArgSort(begin1, end1);
  Reorder(begin2, indices);
}

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

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
    flush();
  }
  void flush() {
    auto str = Gather(buf_.str(), comm_);
    if (MpiWrapper::IsRoot(comm_)) {
      out_ << str;
      out_.flush();
    }
    buf_.str({});
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

template <class M, class T>
std::vector<std::string> FieldToStrings(
    const FieldCell<T>& fcu, const M& m, std::string fmt = "{:4.1f}") {
  using MIdx = typename M::MIdx;
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

template <class M, class T>
void PrintField(
    std::ostream& out, const FieldCell<T>& fcu, M& m,
    std::string fmt = "{:4.1f}") {
  using MIdx = typename M::MIdx;
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

template <class M>
class CommMap {
 public:
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;

  struct State {
    FieldCell<Scal> fc_rank; // rank owning a cell:
    // inner cell: current rank
    // halo cell:
    // * inside the domain: rank of intermediate neighbor
    // * outside the domain:
    //   * lower: rank of neighbor through periodic boundary
    //            minus the number of ranks
    //   * upper: rank of neighbor through periodic boundary
    //            plus the number of ranks
    FieldCell<Scal> fc_col;
    FieldCell<bool> fc_has_neighbors; // true if cell has remote neighbors
    std::vector<IdxCell> flat_cells;
    std::vector<int> flat_cols;
    std::vector<std::array<IdxNci, M::kCellNumNeighborFaces>> flat_nci;
    // indices in `flat_cells` corresponding to inner cells
    // 0: without remote neighbors,  1: with remote neighbors
    std::array<GRange<size_t>, 2> range_inner;
    // indices in `flat_cells` corresponding to halo cells
    GRange<size_t> range_halo;
  };

  struct System {
    bool islead = false; // true on leading block
    int n; // number of rows
    int nnz; // number of non-zero elements
    std::vector<Scal> data;
    std::vector<int> cols;
    std::vector<int> row_ptrs;
    std::map<int, std::vector<int>> recv;
    std::map<int, std::vector<int>> send;

    std::vector<int> neighbors;
    std::vector<int> send_sizes;
    std::vector<const int*> send_maps;
    std::vector<int> recv_sizes;
    std::vector<const int*> recv_maps;

    std::vector<Scal> rhs_data;
    std::vector<Scal> sol_data;
  };
  struct Buffers {
    Scal* rhs;
    Scal* sol;
    size_t size;
  };

  // Returns current System instance valid on leading block.
  // Initialized by Init() and ConvertSystem()
  const System& GetSystem() const {
    fassert(system_.islead, "GetSystem() is only valid on a leading block");
    return system_;
  }
  // Returns pointers to buffers for RHS and solution
  // shared by all local blocks.
  const Buffers& GetBuffers() const {
    fassert(init_called_, "Call Init() before calling GetBuffers()");
    return buffers_;
  }
  void Init(M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      int col_current = 0;
      // Pointers to objects from other local blocks:
      std::vector<const M*> meshes;
      std::vector<State*> states;
    } * ctx(sem);
    auto& t = *ctx;

    const auto comm = m.GetMpiComm();
    const int rank = MpiWrapper::GetCommRank(comm);
    const int commsize = MpiWrapper::GetCommSize(comm);
    if (sem()) {
      auto& s = state_;
      t.meshes.push_back(&m);
      t.states.push_back(&s);
      m.GatherToLead(&t.meshes);
      m.GatherToLead(&t.states);

      s.fc_col.Reinit(m, -1);
      s.fc_rank.Reinit(m, 0);
      for (auto c : m.Cells()) {
        s.fc_rank[c] = rank;
      }
      m.Comm(&s.fc_rank);
    }
    if (sem()) {
      auto& s = state_;
      // Mark ranks through periodic boundaries
      const auto globalsize = m.GetGlobalSize();
      for (auto c : m.SuCellsM()) {
        for (auto d : GRange<size_t>(M::dim)) {
          if (s.fc_rank[c] != rank) {
            if (c[d] < 0) {
              s.fc_rank[c] -= commsize;
            }
            if (c[d] >= globalsize[d]) {
              s.fc_rank[c] += commsize;
            }
          }
        }
      }

      s.fc_has_neighbors.Reinit(m, false);
      for (auto c : m.Cells()) {
        for (auto q : m.Nci(c)) {
          const auto cn = m.GetCell(c, q);
          if (s.fc_rank[cn] != rank) {
            s.fc_has_neighbors[c] = true;
          }
        }
      }
    }
    if (sem("col-inner-no-neghbors")) {
      if (m.IsLead()) {
        t.col_current = 0;
        for (auto part : {0, 1}) {
          for (size_t b = 0; b < t.meshes.size(); ++b) {
            auto& mb = *t.meshes[b];
            auto& sb = *t.states[b];
            const size_t begin = sb.flat_cells.size();
            for (auto c : mb.Cells()) {
              // cells that do not have neighbors from remote ranks
              if (sb.fc_has_neighbors[c] == part) {
                sb.fc_col[c] = t.col_current;
                sb.flat_cells.push_back(c);
                sb.flat_cols.push_back(t.col_current);
                sb.flat_nci.emplace_back();
                ++t.col_current;
                auto& nci = sb.flat_nci.back();
                size_t j = 0;
                for (auto q : mb.Nci(c)) {
                  const auto cn = mb.GetCell(c, q);
                  if (!sb.fc_has_neighbors[cn]) {
                    nci[j++] = q;
                  }
                }
                for (auto q : mb.Nci(c)) {
                  const auto cn = mb.GetCell(c, q);
                  if (sb.fc_has_neighbors[cn]) {
                    nci[j++] = q;
                  }
                }
              }
            }
            sb.range_inner[part] = {begin, sb.flat_cells.size()};
          }
        }
      }
    }
    if (sem("col-inner-other")) {
      // fill halo cells with already computed indices
      m.Comm(&state_.fc_col);
    }
    if (sem("col-halo")) {
      if (m.IsLead()) {
        for (size_t b = 0; b < t.meshes.size(); ++b) {
          auto& mb = *t.meshes[b];
          auto& sb = *t.states[b];
          for (auto c : mb.Cells()) {
            // halo cells from remote ranks
            // -1 to prevent traversing the same cell twice
            for (auto q : mb.Nci(c)) {
              const auto cn = mb.GetCell(c, q);
              if (sb.fc_rank[cn] != rank) {
                sb.fc_col[cn] = -1;
              }
            }
          }
          const size_t begin = sb.flat_cells.size();
          for (auto c : mb.Cells()) {
            // halo cells from remote ranks
            for (auto q : mb.Nci(c)) {
              const auto cn = mb.GetCell(c, q);
              if (sb.fc_rank[cn] != rank) {
                if (sb.fc_col[cn] == -1) {
                  sb.fc_col[cn] = t.col_current;
                  sb.flat_cells.push_back(cn);
                  sb.flat_cols.push_back(t.col_current);
                  ++t.col_current;
                }
              }
            }
          }
          sb.range_halo = {begin, sb.flat_cells.size()};
        }
      }
    }
    if (sem("post") && m.IsLead()) {
      system_.islead = true;
      for (size_t b = 0; b < t.meshes.size(); ++b) {
        auto& mb = *t.meshes[b];
        auto& sb = *t.states[b];
        // collect halo cells to receive
        for (size_t i : sb.range_halo) {
          const auto c = sb.flat_cells[i];
          system_.recv[sb.fc_rank[c]].push_back(sb.flat_cols[i]);
        }
        // collect inner cells to send
        for (size_t i : sb.range_inner[1]) {
          const auto c = sb.flat_cells[i];
          for (auto q : mb.Nci(c)) {
            const auto cn = mb.GetCell(c, q);
            if (sb.fc_rank[cn] != rank) {
              system_.send[sb.fc_rank[cn]].push_back(sb.flat_cols[i]);
            }
          }
        }
      }

      auto realrank = [commsize](int r) {
        while (r < 0) {
          r += commsize;
        }
        while (r >= commsize) {
          r -= commsize;
        }
        return r;
      };

      // merge indices of inner cells to send
      // changing the order of neighbors through periodic boundaries
      for (auto it = system_.send.begin(); it != system_.send.end();) {
        const int r = it->first;
        auto& vreal = system_.send[realrank(r)];
        auto& v = system_.send[r];
        if (r < 0) {
          vreal.insert(vreal.begin(), v.begin(), v.end());
          it = system_.send.erase(it);
        } else if (r >= commsize) {
          vreal.insert(vreal.end(), v.begin(), v.end());
          it = system_.send.erase(it);
        } else {
          ++it;
        }
      }

      // merge indices of halo cells to receive
      // changing the order of neighbors through periodic boundaries
      for (auto it = system_.recv.begin(); it != system_.recv.end();) {
        const int r = it->first;
        auto& vreal = system_.recv[realrank(r)];
        auto& v = system_.recv[r];
        if (r < 0) {
          vreal.insert(vreal.end(), v.begin(), v.end());
          it = system_.recv.erase(it);
        } else if (r >= commsize) {
          vreal.insert(vreal.begin(), v.begin(), v.end());
          it = system_.recv.erase(it);
        } else {
          ++it;
        }
      }

      std::vector<int> send_neighbors;
      for (const auto& p : system_.send) {
        send_neighbors.push_back(p.first);
        system_.send_sizes.push_back(p.second.size());
        system_.send_maps.push_back(p.second.data());
      }

      std::vector<int> recv_neighbors;
      for (const auto& p : system_.recv) {
        recv_neighbors.push_back(p.first);
        system_.recv_sizes.push_back(p.second.size());
        system_.recv_maps.push_back(p.second.data());
      }

      fassert_equal(send_neighbors, recv_neighbors);
      system_.neighbors = send_neighbors;
    }
    if (sem()) {
      init_called_ = true;
    }
  }
  void ConvertSystem(const FieldCell<Expr>& fc_system, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      // Pointers to objects from other local blocks:
      std::vector<const M*> meshes;
      std::vector<State*> states;
      std::vector<const FieldCell<Expr>*> fields;
      std::vector<Buffers*> buffers;
      int col_current = 0;
    } * ctx(sem);
    auto& t = *ctx;

    if (sem("gather")) {
      fassert(init_called_, "Call Init() before calling ConvertSystem()");
      t.meshes.push_back(&m);
      t.states.push_back(&state_);
      t.fields.push_back(&fc_system);
      t.buffers.push_back(&buffers_);
      m.GatherToLead(&t.meshes);
      m.GatherToLead(&t.states);
      m.GatherToLead(&t.fields);
      m.GatherToLead(&t.buffers);
    }
    if (sem("data") && m.IsLead()) {
      system_.cols.clear();
      system_.data.clear();
      system_.row_ptrs = {0};

      for (auto part : {0, 1}) {
        for (size_t b = 0; b < t.meshes.size(); ++b) {
          auto& mb = *t.meshes[b];
          auto& sb = *t.states[b];
          auto& fb = *t.fields[b];
          for (size_t i : sb.range_inner[part]) {
            const auto c = sb.flat_cells[i];
            const int col_diag = sb.fc_col[c];
            Scal data_diag = fb[c][0];
            for (auto q : sb.flat_nci[i]) {
              const auto cn = mb.GetCell(c, q);
              const int col = sb.fc_col[cn];
              if (col == col_diag) { // reduce repeting diagonal terms (for 2D)
                data_diag += fb[c][1 + q.raw()];
              }
            }
            system_.cols.push_back(sb.fc_col[c]);
            system_.data.push_back(data_diag);
            for (auto q : sb.flat_nci[i]) {
              const auto cn = mb.GetCell(c, q);
              const int col = sb.fc_col[cn];
              if (col != col_diag) {
                system_.cols.push_back(col);
                system_.data.push_back(fb[c][1 + q.raw()]);
              }
            }
            system_.row_ptrs.push_back(system_.cols.size());
          }
        }
      }

      fassert(system_.row_ptrs.size() >= 1);
      system_.n = system_.row_ptrs.size() - 1;
      system_.nnz = system_.data.size();

      // allocate buffers and copy pointers to all blocks
      system_.rhs_data.resize(system_.n);
      system_.sol_data.resize(system_.n);
      for (size_t b = 0; b < t.meshes.size(); ++b) {
        auto& bb = *t.buffers[b];
        bb.rhs = system_.rhs_data.data();
        bb.sol = system_.sol_data.data();
        bb.size = system_.n;
      }
    }
  }
  void PrintStat(M& m, std::ostream& outstream) const {
    auto sem = m.GetSem(__func__);
    struct {
      FieldCell<int> fc_nci;
    } * ctx(sem);
    const auto comm = m.GetMpiComm();
    const int rank = MpiWrapper::GetCommRank(comm);
    auto& t = *ctx;
    auto& s = state_;
    if (sem()) {
      t.fc_nci.Reinit(m, 0);
      for (size_t i = 0; i < s.flat_nci.size(); ++i) {
        const auto c = s.flat_cells[i];
        Scal w = 0;
        for (auto q : s.flat_nci[i]) {
          w = w * 10 + q.raw();
        }
        t.fc_nci[c] = w;
      }
    }
    if (sem.Nested()) {
      PrintField(outstream, s.fc_rank, m, "{:2g}");
    }
    if (sem.Nested()) {
      PrintField(outstream, s.fc_col, m, "{:3g}");
    }
    if (sem.Nested()) {
      PrintField(outstream, t.fc_nci, m, " {:06d}");
    }
    if (sem() && m.IsLead()) {
      auto& system = GetSystem();
      StreamMpi out(outstream, comm);
      out << "\nrank=" << rank;
      for (const auto& p : system.send) {
        out << "\nsend to rank " << p.first;
        out << "\n" << p.second;
      }
      for (const auto& p : system.recv) {
        out << "\nrecv from rank " << p.first;
        out << "\n" << p.second;
      }
      out << "\nsystem.data";
      out << "\n" << system.data;

      out << "\nsystem.cols";
      out << "\n" << system.cols;

      out << "\nsystem.row_ptrs";
      out << "\n" << system.row_ptrs;

      out << "\nsystem.neighbors";
      out << "\n" << system.neighbors;

      out << "\nsystem.send_sizes";
      out << "\n" << system.send_sizes;

      out << "\nsystem.recv_sizes";
      out << "\n" << system.recv_sizes;
    }
  }
  void PrintStat(M& m) const {
    PrintStat(m, std::cerr);
  }
  // Copies field to a flat array.
  void FieldToArray(const FieldCell<Scal>& fcu, Scal* buf) const {
    fassert(buf);
    const auto& s = state_;
    for (auto range : s.range_inner) {
      for (size_t i : range) {
        buf[s.flat_cols[i]] = fcu[s.flat_cells[i]];
      }
    }
  }
  // Copied flat array to fields.
  // flat: pointer to array, significant only on leading block.
  // Returns field with filled inner cells.
  void ArrayToField(const Scal* buf, FieldCell<Scal>& fcu, const M& m) const {
    const auto& s = state_;
    fcu.Reinit(m);
    for (auto range : s.range_inner) {
      for (size_t i : range) {
        fcu[s.flat_cells[i]] = buf[s.flat_cols[i]];
      }
    }
  }
  void ArrayToField(
      const std::vector<Scal>& buf, FieldCell<Scal>& fcu, const M& m) const {
    ArrayToField(buf.data(), fcu, m);
  }

 private:
  System system_; // valid on leading block, empty on others
  State state_;
  Buffers buffers_;
  bool init_called_ = false;
};
