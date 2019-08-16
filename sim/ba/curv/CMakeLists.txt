cmake_minimum_required(VERSION 2.8)
set(CHPREFIX $ENV{CHPREFIX})
project(curv NONE)

set(T curv)
set(D ${CHPREFIX}/include)
set(H
  curv/div.h
  curv/fix.h
  curv/learn.h
  curv/partstr.h
  curv/select.h
  curv/zero.h
  curv/learn_eval.h
)

configure_file(conf.in ${T}.conf @ONLY)
set(DD ${CHPREFIX}/include/curv)
install(FILES ${H} DESTINATION ${DD})
install(FILES  ${T}.conf DESTINATION ${DD})
