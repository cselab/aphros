set(T inside)
# FIXME: these files should be compiled to a single object library
#        but CMake creates one object for every source file
#        so they a compiled to a static library instead
#        and linked to aphros through a separate list BinLibs
add_library(${T} STATIC
    main.c
    bbox.c
    err.c
    memory.c
    off.c
    ply.c
    predicate.c
    predicate/main.c
    stl.c
  )
