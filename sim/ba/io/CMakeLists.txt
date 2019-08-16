cmake_minimum_required(VERSION 2.8)
set(CHPREFIX $ENV{CHPREFIX})
project(io NONE)

set(T io)
set(D ${CHPREFIX}/include)
set(H
  io/2d.h
  io/3d.h
  io/io.h
  io/iompi.h
  io/util.h
  io/write.h
)

configure_file(conf.in ${T}.conf @ONLY)
set(DD ${CHPREFIX}/include/io)
install(FILES ${H} DESTINATION ${DD})
install(FILES  ${T}.conf DESTINATION ${DD})
