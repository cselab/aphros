#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

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

template <class M>
class CommMap {
 public:
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  void Run(M& m, const FieldCell<Expr>& fc_system) {
    auto comm = m.GetMpiComm();

    FieldCell<int> cell_to_rank(m);
    int rank;
    MPI_Comm_rank(comm, &rank);
    for (auto c : m.Cells()) {
      cell_to_rank[c] = rank;
    }
    m.Comm(&cell_to_rank);

    std::vector<IdxCell> row_to_cell;
    std::vector<char> row_to_num_neighbors;
    for (auto c : m.Cells()) {
      row_to_cell.push_back(c);
      row_to_num_neighbors.push_back(0);
      for (auto q : m.Nci(c)) {
        if (cell_to_rank[m.GetCell(c, q)] != rank) {
          ++row_to_num_neighbors.back();
        }
      }
    }

    {
      const auto indices =
          ArgSort(row_to_num_neighbors.begin(), row_to_num_neighbors.end());
      Reorder(row_to_cell.begin(), indices);
      Reorder(row_to_num_neighbors.begin(), indices);
    }

    int commsize;
    MPI_Comm_size(comm, &commsize);

    FieldCell<int> cell_to_rankext = cell_to_rank;
    for (auto c : m.AllCells()) {
      if (cell_to_rankext[c] != rank) {
        cell_to_rankext[c] += commsize;
      }
    }

    std::vector<IdxCell> flat_to_cell;
    std::vector<char> flat_to_coeff;
    std::vector<int> flat_to_rankext;
    for (auto c : row_to_cell) {
      flat_to_cell.push_back(c);
      flat_to_coeff.push_back(0);
      flat_to_rankext.push_back(cell_to_rankext[c]);
      for (auto q : m.Nci(c)) {
        flat_to_cell.push_back(c);
        flat_to_coeff.push_back(q + 1);
        flat_to_rankext.push_back(cell_to_rankext[m.GetCell(c, q)]);
      }
    }

    {
      const auto indices =
          ArgSort(flat_to_rankext.begin(), flat_to_rankext.end());
      Reorder(flat_to_cell.begin(), indices);
      Reorder(flat_to_rankext.begin(), indices);
      Reorder(flat_to_coeff.begin(), indices);
    }

    std::vector<Scal> data(flat_to_cell.size());
    std::vector<int> row_ptrs(flat_to_cell.size());

    FieldCell<int> cell_to_col(m, -1);
    {
      int col = 0;
      for (auto c : row_to_cell) {
        cell_to_col[c] = col++;
      }
      for (size_t i = 0; i < flat_to_cell.size(); ++i) {
        const int q = flat_to_coeff[i] - 1;
        if (q >= 0) {
          const auto c = m.GetCell(flat_to_cell[i], q);
          if (cell_to_col[c] < 0) {
            cell_to_col[c] = col++;
          }
        }
      }
    }

    std::vector<int> recv_ptrs;
    std::vector<int> recv_maps;
    IdxCell cprev = IdxCell(-1);
    for (size_t i = 0; i < flat_to_cell.size(); ++i) {
      const auto c = flat_to_cell[i];
      data[i] = fc_system[c][flat_to_coeff[i]];
      if (c != cprev) {
        row_ptrs.push_back(cell_to_col.size());
      }
      cprev = c;
    }
  }
};

