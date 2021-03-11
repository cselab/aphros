#include "system.h"
#ifdef _WIN32
#include "system_windows.ipp"
#else
#include "system_unix.ipp"
#endif
