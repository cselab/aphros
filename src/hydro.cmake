include("$ENV{CHPREFIX}/include/client.cmake")

set(T "hydro")
add_library(${T} INTERFACE)
link_directories(${CHPREFIX}/lib)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_LINK_LIBRARIES distr local cubism hypre report suspender sysinfo utilfluid events timer git
gitgen simple vof normal partstrmesh tvd convdiffvi convdiffi convdiffe solver
parser vars dumper overlap young hypreext hdf)
