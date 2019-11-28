#include "subcomm.h"
#include "sysinfo.h"

// Example for two nodes each with two threads.
// Ranks:
// comm_world:  | 0 1 | 2 3 |
// comm_omp:    | 0 1 | 0 1 |
// comm_master: | 0 - | 1 - |



