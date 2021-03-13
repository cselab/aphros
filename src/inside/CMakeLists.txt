set(T inside)
# FIXME: these files should be compiled to a single object library
#        but CMake creates one object for every source file
#        so here they are compiled to a static library instead
#        and linked to aphros through a separate list StaticLibs
add_library(${T} STATIC
    main.c
    bbox.c
    err.c
    memory.c
    off.c
    ply.c
    predicate.c
    stl.c
  )
install(TARGETS ${T} EXPORT export DESTINATION lib)
