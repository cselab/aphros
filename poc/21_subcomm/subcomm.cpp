#include <algorithm>
#include <cassert>
#include <map>
#include <mpi.h>
#include <set>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm comm_hypre_ = MPI_COMM_WORLD;
    MPI_Comm comm_, comm_omp_;
    MPI_Group omp_master_group_;
    bool omp_master_;

    ////////////////////////////////////////////////////////////////////////////
    // 1. Find color of ranks in comm_hypre_ that belong to same socket and
    // create split communicator stored in comm_omp_
    // 2. Create new communicator with only root ranks in comm_omp_, store in
    // temporarily in comm_

    // 1.
    unsigned long a, d, c;
    __asm__ volatile("rdtscp" : "=a"(a), "=d"(d), "=c"(c));
    const size_t node_ID = (c & 0xFFF000) >> 12;

    std::string hostname(1024, '\0'); // max. hostname length (Linux unistd.h)
    gethostname(&hostname.front(), hostname.size());
    const size_t host_ID = std::hash<std::string>{}(hostname);
    const size_t c_hash = host_ID ^ (node_ID << 1);

    // construct color and split communicator
    int hypre_size, hypre_rank;
    MPI_Comm_size(comm_hypre_, &hypre_size);
    MPI_Comm_rank(comm_hypre_, &hypre_rank);
    std::vector<size_t> all_hashes(hypre_size);
    MPI_Allgather(&c_hash,
                  1,
                  MPI_UINT64_T,
                  all_hashes.data(),
                  1,
                  MPI_UINT64_T,
                  comm_hypre_);
    std::set<size_t> unique_nodes(all_hashes.begin(), all_hashes.end());
    std::map<size_t, int> cmap;
    int color = 0;
    for (auto ID : unique_nodes) {
        cmap[ID] = color++;
    }
    color = cmap[c_hash];
    MPI_Comm_split(comm_hypre_, color, hypre_rank, &comm_omp_);

    // 2.
    int omp_rank;
    MPI_Comm_rank(comm_omp_, &omp_rank);
    if (0 != omp_rank) {
        omp_rank = -omp_rank; // if I am not root in comm_omp_
        omp_master_ = false;
    } else {
        omp_master_ = true;
    }
    std::vector<int> omp2hypre(hypre_size);
    MPI_Allgather(
        &omp_rank, 1, MPI_INT, omp2hypre.data(), 1, MPI_INT, comm_hypre_);
    std::set<int> m2h(omp2hypre.begin(), omp2hypre.end());
    omp2hypre.clear();
    for (auto r : m2h) {
        if (r >= 0) {
            omp2hypre.push_back(hypre_rank);
        }
    }
    std::sort(omp2hypre.begin(), omp2hypre.end());

    MPI_Group hypre_group;
    MPI_Comm_group(comm_hypre_, &hypre_group);

    MPI_Group_incl(hypre_group,
                   static_cast<int>(omp2hypre.size()),
                   omp2hypre.data(),
                   &omp_master_group_);
    MPI_Comm_create(comm_hypre_, omp_master_group_, &comm_);
    MPI_Group_free(&hypre_group);

    {
        const size_t core_ID = c & 0xFFF;
        int omp_size, comm_rank, comm_size;
        MPI_Comm_size(comm_omp_, &omp_size);
        comm_rank = -1;
        comm_size = -1;
        if (omp_master_) {
            MPI_Comm_size(comm_, &comm_size);
            MPI_Comm_rank(comm_, &comm_rank);
        }
        printf(
            "HYBRID: host:%s; comm_hypre_:%d/%d; comm_omp_:%d/%d; comm_:%d/%d; "
            "node_ID=%lu; core_ID=%lu\n",
            hostname.c_str(),
            hypre_rank,
            hypre_size,
            omp_rank,
            omp_size,
            comm_rank,
            comm_size,
            node_ID,
            core_ID);
    }

    ////////////////////////////////////////////////////////////////////////////
    // do some work
    int rank = -1;
    if (omp_master_) { // only subcomm ranks must do this
        MPI_Comm_rank(comm_, &rank);
    }
    printf("My comm_ rank: %d\n", rank);

    MPI_Finalize();
    return 0;
}
