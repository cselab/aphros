include("$ENV{CHPREFIX}/include/client.cmake")

set(T "hydro")
add_library(${T} INTERFACE)
target_link_directories(${T} PUBLIC INTERFACE ${CHPREFIX}/lib)
target_link_libraries(${T} INTERFACE
distr local cubism hypre report suspender sysinfo utilfluid events timer git
gitgen simple vof normal partstrmesh tvd convdiffvi convdiffi convdiffe solver
parser vars dumper overlap young hypreext hdf)
