set(T inside)
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
