cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphrosjs)

option(BUILD_DIFFUSION "Build difusion demo diffusion.js" ON)
option(BUILD_ELECTROCHEM "Build reactor demo electrochem.js" ON)
option(BUILD_EXPLORER "Build Aphros explorer aphros.js" ON)
option(BUILD_HYDRO "Build Multi-VOF demo hydro.js" ON)

add_link_options(-sALLOW_MEMORY_GROWTH=1)
add_compile_definitions(_USE_MPI_=0)
include_directories(../../src)

if (BUILD_DIFFUSION)
  set(T diffusion)
  add_executable(${T} diffusion.cpp)
  target_link_libraries(${T} aphros)
  target_link_options(${T} PRIVATE "-sEXPORTED_FUNCTIONS=[\
    '_main', '_MulDiffusion', '_AddVelocityAngle', \
    '_TogglePause', '_Init', '_SetMesh' \
  ]")
  target_link_options(${T} PRIVATE "-sEXPORTED_RUNTIME_METHODS=['cwrap', 'ccall']")
  configure_file(diffusion.html diffusion.html COPYONLY)
endif()

if (BUILD_ELECTROCHEM)
  set(T electrochem)
  add_executable(${T} electrochem.cpp)
  target_link_libraries(${T} aphros)
  target_link_options(${T} PRIVATE "-sEXPORTED_FUNCTIONS=[\
    '_main', '_SetRuntimeConfig', '_SetExtraConfig', \
    '_TogglePause', '_Init', '_Spawn', '_SetMesh' ,'_GetLines' \
  ]")
  target_link_options(${T} PRIVATE "-sEXPORTED_RUNTIME_METHODS=['cwrap', 'ccall']")
  configure_file(electrochem.html electrochem.html COPYONLY)
  configure_file(electrochem_defs.js electrochem_defs.js COPYONLY)
endif()

if (BUILD_EXPLORER)
  set(T wexplorer)
  add_executable(${T} aphros.cpp)
  target_link_libraries(${T} aphros)
  target_link_options(${T} PRIVATE --preload-file conf)
  target_link_options(${T} PRIVATE "-sEXPORTED_FUNCTIONS=[\
    '_main', '_GetConfigString', '_GetConfigDouble',\
    '_SetRuntimeConfig', '_SetExtraConfig', '_TogglePause',\
    '_Spawn', '_SetMesh', '_GetLines'\
  ]")
  target_link_options(${T} PRIVATE "-sEXPORTED_RUNTIME_METHODS=['cwrap', 'ccall']")
  set_target_properties(${T} PROPERTIES OUTPUT_NAME aphros)
  configure_file(conf/std.conf conf/std.conf)
  configure_file(../../deploy/scripts/sim_base.conf conf/base.conf COPYONLY)
  configure_file(aphros.html aphros.html COPYONLY)
  configure_file(aphros_defs.js aphros_defs.js COPYONLY)
  configure_file(libs/lz-string.js libs/lz-string.js COPYONLY)
endif()

if (BUILD_HYDRO)
  set(T whydro)
  add_executable(${T} hydro.cpp)
  target_link_libraries(${T} aphros)
  target_link_options(${T} PRIVATE "-sEXPORTED_FUNCTIONS=[\
    '_main', '_AddVelocityAngle', '_SetRuntimeConfig',\
    '_SetExtraConfig', '_TogglePause', '_Init', \
    '_Spawn', '_SetMesh' ,'_GetLines', '_SetCoal' \
  ]")
  target_link_options(${T} PRIVATE "-sEXPORTED_RUNTIME_METHODS=['cwrap', 'ccall']")
  set_target_properties(${T} PROPERTIES OUTPUT_NAME hydro)
  configure_file(hydro.html hydro.html COPYONLY)
  configure_file(hydro_defs.js hydro_defs.js COPYONLY)
endif()

