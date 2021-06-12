B = build

$B/libhook.so: $B/Makefile hook.cpp
	make -C $B

$B/Makefile: CMakeLists.txt
	mkdir -p $B
	cd $B && cmake ..
