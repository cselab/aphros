add_test_current(COMMAND test
  ASSETS
  ${PROJECT_SOURCE_DIR}/../examples/208_rising/plot.py
  ${PROJECT_SOURCE_DIR}/../examples/208_rising/std.conf
  add.conf
  test
  ref
  )
# FIXME: revise without `execute_process`
#        Need to to rename Makefile to sim.makefile
#        `file(COPY_FILE)` not supported in CMake 3.10.2)
execute_process(COMMAND cp -f
  "${PROJECT_SOURCE_DIR}/../examples/208_rising/Makefile"
  "${CMAKE_CURRENT_BINARY_DIR}/sim.makefile")
execute_process(COMMAND mkdir -p "${CMAKE_CURRENT_BINARY_DIR}/ref")
execute_process(COMMAND cp -f
  "${PROJECT_SOURCE_DIR}/../examples/208_rising/ref/c1g3l4s.txt"
  "${CMAKE_CURRENT_BINARY_DIR}/ref/")
