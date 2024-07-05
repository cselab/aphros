set(T "overlap")
add_object(${T} overlap.cpp)
object_link_libraries(${T} use_avx)
