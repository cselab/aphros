BUILD = build
CMAKE = cmake


$(BUILD)/libhook.so: $(BUILD)/Makefile hook.cpp
	make -C $(BUILD)

$(BUILD)/Makefile: CMakeLists.txt
	mkdir -p $(BUILD)
	cd $(BUILD) && $(CMAKE) ..
