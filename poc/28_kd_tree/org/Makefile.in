obj = kdtree.o
name = kdtree
lib_a = lib$(name).a

so_major = 0
so_minor = 1

CFLAGS = -pedantic -Wall $(pic) $(opt) $(dbg) $(falloc) $(pthreads)
LDFLAGS = $(ldpthread)

sys = $(shell uname -s | sed 's/MING32.*/MINGW32/')
ifeq ($(sys), Darwin)
	lib_so = lib$(name).dylib
	shared = -dynamiclib
else ifeq ($(sys), MINGW32)
	lib_so = lib$(name).dll
	shared = -shared -Wl,--output-def,lib$(name).def
else
	devname = lib$(name).so
	soname = lib$(name).so.$(so_major)
	lib_so = lib$(name).so.$(so_major).$(so_minor)
	shared = -shared -Wl,-soname=$(soname)
	pic = -fPIC
endif

.PHONY: all
all: $(lib_a) $(lib_so)

$(lib_a): $(obj)
	$(AR) rcs $@ $(obj)

$(lib_so): $(obj)
	$(CC) $(shared) -o $@ $(obj) $(LDFLAGS)
	[ -n "$(soname)" ] && ln -s $(lib_so) $(soname) || true
	[ -n "$(devname)" ] && ln -s $(soname) $(devname) || true

.PHONY: examples
examples:
	cd examples; $(MAKE)

.PHONY: clean
clean:
	rm -f $(obj) $(lib_a) $(lib_so)

.PHONY: install
install:
	mkdir -p $(DESTDIR)$(PREFIX)/include $(DESTDIR)$(PREFIX)/lib
	cp kdtree.h $(DESTDIR)$(PREFIX)/include/kdtree.h
	cp -d $(lib_so) $(soname) $(devname) $(DESTDIR)$(PREFIX)/lib/
	cp $(lib_a) $(DESTDIR)$(PREFIX)/lib/

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/include/kdtree.h
	rm -f $(DESTDIR)$(PREFIX)/lib/$(lib_so)
	rm -f $(DESTDIR)$(PREFIX)/lib/$(lib_a)
	[ -n "$(soname)" ] && rm -f $(DESTDIR)$(PREFIX)/lib/$(soname) || true
	[ -n "$(devname)" ] && rm -f $(DESTDIR)$(PREFIX)/lib/$(devname) || true
