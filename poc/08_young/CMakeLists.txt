cmake_minimum_required(VERSION 3.12)
project(young, C)
add_library(young STATIC main.c)
set_target_properties(young PROPERTIES PUBLIC_HEADER young.h)
install(TARGETS young
	ARCHIVE DESTINATION lib
	PUBLIC_HEADER DESTINATION include)
