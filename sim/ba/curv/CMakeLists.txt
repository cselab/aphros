cmake_minimum_required(VERSION 2.8)
set(CHPREFIX $ENV{CHPREFIX})
project(curv NONE)

set(T curv)
set(D ${CHPREFIX}/include)
set(H
  curv/curvature_div.h
  curv/curvature_fix.h
  curv/curvature_learn.h
  curv/curvature_partstr.h
  curv/curvature_select.h
  curv/curvature_zero.h
  curv/learn.h
)

configure_file(conf.in ${T}.conf @ONLY)
set(DD ${CHPREFIX}/include/curv)
install(FILES ${H} DESTINATION ${DD})
install(FILES  ${T}.conf DESTINATION ${DD})
