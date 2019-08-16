cmake_minimum_required(VERSION 2.8)
set(CHPREFIX $ENV{CHPREFIX})
project(curv NONE)

set(T curv)
set(D ${CHPREFIX}/include/curv)
set(H
  curvature_div.h
  curvature_fix.h
  curvature_learn.h
  curvature_partstr.h
  curvature_select.h
  curvature_zero.h
  learn.h
)

configure_file(conf.in ${T}.conf @ONLY)
install(FILES ${H} DESTINATION ${D})
install(FILES  ${T}.conf DESTINATION ${D})
