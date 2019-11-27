#include <cstddef>

namespace sysinfo {

// Returns resident memory [bytes]
size_t GetMem();
// Returns true if Intel Hyperthreads are enabled
bool HasHyperthreads();

} // namespace sysinfo
