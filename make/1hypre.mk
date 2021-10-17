CXXFLAGS_HYPRE = -I$(APHROS_PREFIX)/include
LDFLAGS_HYPRE = $(APHROS_PREFIX)/lib/libHYPRE.so
O_HYPRE = \
$(WRK)/linear/hypre.o \
$(WRK)/linear/hypresub.o \
$(WRK)/linear/linear_hypre.o \

