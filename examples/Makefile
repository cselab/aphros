# Examples that only need to be compiled.
DIRS_BUILD = \
	000_march \
	001_march_cell \
	002_march_sphere \
	006_cmake \
	100_advection \
	101_partstr \
	102_mean_curv_flow \
	103_embed_advection \
	104_embed_diffusion \
	105_sharp \
	106_diffusion_4d \

# Examples that run a simulation to build.
DIRS_SIM = \
	200_primlist \
	201_taylor_couette \
	203_hydrostatic \
	204_microfoam \

DIRS = $(DIRS_BUILD) $(DIRS_SIM)

VERBOSE =
REDIR = $(if $(VERBOSE),,>/dev/null 2>&1)
CMAKE = cmake

all:
	+make build
	+make run

build: $(DIRS)

build_nosim: $(DIRS_BUILD)

$(DIRS):
	@echo "$@"
	@(cd "$@" && test -f Makefile || ! test -f CMakeLists.txt ||  $(CMAKE) . ) $(REDIR)
	@+make -s -C "$@" $(REDIR)

run:
	@for d in $(DIRS) ; do echo "$$d" ; (cd "$$d" && ./run $(REDIR) ) ; done

clean:
	@for d in $(DIRS) ; do echo "$$d" ; make -C "$$d" clean $(REDIR) ; done

.PHONY: all build build_nosim run clean $(DIRS)
