$(WRK)/aphros_c/git.o: $(SRC)/aphros_c/git.cpp; $(CXX_RULE) $(SRC)/aphros_c/git.cpp
$(WRK)/aphros_c/main.o: $(SRC)/aphros_c/main.cpp; $(CXX_RULE) $(SRC)/aphros_c/main.cpp
$(WRK)/aphros_c/parser.o: $(SRC)/aphros_c/parser.cpp; $(CXX_RULE) $(SRC)/aphros_c/parser.cpp
$(WRK)/color/color.o: $(SRC)/color/color.c; $(CC_RULE) $(SRC)/color/color.c
$(WRK)/distr/comm_manager.o: $(SRC)/distr/comm_manager.cpp; $(CXX_RULE) $(SRC)/distr/comm_manager.cpp
$(WRK)/distr/distrbasic.o: $(SRC)/distr/distrbasic.cpp; $(CXX_RULE) $(SRC)/distr/distrbasic.cpp
$(WRK)/distr/distr.o: $(SRC)/distr/distr.cpp; $(CXX_RULE) $(SRC)/distr/distr.cpp
$(WRK)/distr/distr_particles.o: $(SRC)/distr/distr_particles.cpp; $(CXX_RULE) $(SRC)/distr/distr_particles.cpp
$(WRK)/distr/distrsolver.o: $(SRC)/distr/distrsolver.cpp; $(CXX_RULE) $(SRC)/distr/distrsolver.cpp
$(WRK)/distr/local.o: $(SRC)/distr/local.cpp; $(CXX_RULE) $(SRC)/distr/local.cpp
$(WRK)/distr/native.o: $(SRC)/distr/native.cpp; $(CXX_RULE) $(SRC)/distr/native.cpp
$(WRK)/distr/report.o: $(SRC)/distr/report.cpp; $(CXX_RULE) $(SRC)/distr/report.cpp
$(WRK)/dump/dump.o: $(SRC)/dump/dump.cpp; $(CXX_RULE) $(SRC)/dump/dump.cpp
$(WRK)/dump/dumper.o: $(SRC)/dump/dumper.cpp; $(CXX_RULE) $(SRC)/dump/dumper.cpp
$(WRK)/dump/hdf.o: $(SRC)/dump/hdf.cpp; $(CXX_RULE) $(SRC)/dump/hdf.cpp
$(WRK)/dump/raw.o: $(SRC)/dump/raw.cpp; $(CXX_RULE) $(SRC)/dump/raw.cpp
$(WRK)/dump/xmf.o: $(SRC)/dump/xmf.cpp; $(CXX_RULE) $(SRC)/dump/xmf.cpp
$(WRK)/explorer.o: $(SRC)/explorer.cpp; $(CXX_RULE) $(SRC)/explorer.cpp
$(WRK)/func/init_contang.o: $(SRC)/func/init_contang.cpp; $(CXX_RULE) $(SRC)/func/init_contang.cpp
$(WRK)/func/init.o: $(SRC)/func/init.cpp; $(CXX_RULE) $(SRC)/func/init.cpp
$(WRK)/func/init_vel.o: $(SRC)/func/init_vel.cpp; $(CXX_RULE) $(SRC)/func/init_vel.cpp
$(WRK)/func/primlist.o: $(SRC)/func/primlist.cpp; $(CXX_RULE) $(SRC)/func/primlist.cpp
$(WRK)/geom/mesh.o: $(SRC)/geom/mesh.cpp; $(CXX_RULE) $(SRC)/geom/mesh.cpp
$(WRK)/inside/bbox.o: $(SRC)/inside/bbox.c; $(CC_RULE) $(SRC)/inside/bbox.c
$(WRK)/inside/err.o: $(SRC)/inside/err.c; $(CC_RULE) $(SRC)/inside/err.c
$(WRK)/inside/main.o: $(SRC)/inside/main.c; $(CC_RULE) $(SRC)/inside/main.c
$(WRK)/inside/memory.o: $(SRC)/inside/memory.c; $(CC_RULE) $(SRC)/inside/memory.c
$(WRK)/inside/off.o: $(SRC)/inside/off.c; $(CC_RULE) $(SRC)/inside/off.c
$(WRK)/inside/ply.o: $(SRC)/inside/ply.c; $(CC_RULE) $(SRC)/inside/ply.c
$(WRK)/inside/predicate.o: $(SRC)/inside/predicate.c; $(CC_RULE) $(SRC)/inside/predicate.c
$(WRK)/inside/stl.o: $(SRC)/inside/stl.c; $(CC_RULE) $(SRC)/inside/stl.c
$(WRK)/kernel/hydro.o: $(SRC)/kernel/hydro.cpp; $(CXX_RULE) $(SRC)/kernel/hydro.cpp
$(WRK)/linear/linear.o: $(SRC)/linear/linear.cpp; $(CXX_RULE) $(SRC)/linear/linear.cpp
$(WRK)/main.o: $(SRC)/main.c; $(CC_RULE) $(SRC)/main.c
$(WRK)/march/main.o: $(SRC)/march/main.c; $(CC_RULE) $(SRC)/march/main.c
$(WRK)/overlap/overlap.o: $(SRC)/overlap/overlap.cpp; $(CXX_RULE) $(SRC)/overlap/overlap.cpp
$(WRK)/parse/argparse.o: $(SRC)/parse/argparse.cpp; $(CXX_RULE) $(SRC)/parse/argparse.cpp
$(WRK)/parse/codeblocks.o: $(SRC)/parse/codeblocks.cpp; $(CXX_RULE) $(SRC)/parse/codeblocks.cpp
$(WRK)/parse/conf2py.o: $(SRC)/parse/conf2py.cpp; $(CXX_RULE) $(SRC)/parse/conf2py.cpp
$(WRK)/parse/parser.o: $(SRC)/parse/parser.cpp; $(CXX_RULE) $(SRC)/parse/parser.cpp
$(WRK)/parse/template.o: $(SRC)/parse/template.cpp; $(CXX_RULE) $(SRC)/parse/template.cpp
$(WRK)/parse/vars.o: $(SRC)/parse/vars.cpp; $(CXX_RULE) $(SRC)/parse/vars.cpp
$(WRK)/solver/approx_eb.o: $(SRC)/solver/approx_eb.cpp; $(CXX_RULE) $(SRC)/solver/approx_eb.cpp
$(WRK)/solver/approx.o: $(SRC)/solver/approx.cpp; $(CXX_RULE) $(SRC)/solver/approx.cpp
$(WRK)/solver/convdiffe.o: $(SRC)/solver/convdiffe.cpp; $(CXX_RULE) $(SRC)/solver/convdiffe.cpp
$(WRK)/solver/convdiffi.o: $(SRC)/solver/convdiffi.cpp; $(CXX_RULE) $(SRC)/solver/convdiffi.cpp
$(WRK)/solver/convdiffvg.o: $(SRC)/solver/convdiffvg.cpp; $(CXX_RULE) $(SRC)/solver/convdiffvg.cpp
$(WRK)/solver/curv.o: $(SRC)/solver/curv.cpp; $(CXX_RULE) $(SRC)/solver/curv.cpp
$(WRK)/solver/electro.o: $(SRC)/solver/electro.cpp; $(CXX_RULE) $(SRC)/solver/electro.cpp
$(WRK)/solver/embed.o: $(SRC)/solver/embed.cpp; $(CXX_RULE) $(SRC)/solver/embed.cpp
$(WRK)/solver/fluid_dummy.o: $(SRC)/solver/fluid_dummy.cpp; $(CXX_RULE) $(SRC)/solver/fluid_dummy.cpp
$(WRK)/solver/normal.o: $(SRC)/solver/normal.cpp; $(CXX_RULE) $(SRC)/solver/normal.cpp
$(WRK)/solver/particles.o: $(SRC)/solver/particles.cpp; $(CXX_RULE) $(SRC)/solver/particles.cpp
$(WRK)/solver/partstrmeshm.o: $(SRC)/solver/partstrmeshm.cpp; $(CXX_RULE) $(SRC)/solver/partstrmeshm.cpp
$(WRK)/solver/proj_eb.o: $(SRC)/solver/proj_eb.cpp; $(CXX_RULE) $(SRC)/solver/proj_eb.cpp
$(WRK)/solver/proj.o: $(SRC)/solver/proj.cpp; $(CXX_RULE) $(SRC)/solver/proj.cpp
$(WRK)/solver/simple.o: $(SRC)/solver/simple.cpp; $(CXX_RULE) $(SRC)/solver/simple.cpp
$(WRK)/solver/solver.o: $(SRC)/solver/solver.cpp; $(CXX_RULE) $(SRC)/solver/solver.cpp
$(WRK)/solver/tracer.o: $(SRC)/solver/tracer.cpp; $(CXX_RULE) $(SRC)/solver/tracer.cpp
$(WRK)/solver/vofm.o: $(SRC)/solver/vofm.cpp; $(CXX_RULE) $(SRC)/solver/vofm.cpp
$(WRK)/solver/vof.o: $(SRC)/solver/vof.cpp; $(CXX_RULE) $(SRC)/solver/vof.cpp
$(WRK)/test/advection/main.o: $(SRC)/test/advection/main.cpp; $(CXX_RULE) $(SRC)/test/advection/main.cpp
$(WRK)/test/approx/linear.o: $(SRC)/test/approx/linear.cpp; $(CXX_RULE) $(SRC)/test/approx/linear.cpp
$(WRK)/test/approx/main.o: $(SRC)/test/approx/main.cpp; $(CXX_RULE) $(SRC)/test/approx/main.cpp
$(WRK)/test/benchmark_diffusion/main.o: $(SRC)/test/benchmark_diffusion/main.cpp; $(CXX_RULE) $(SRC)/test/benchmark_diffusion/main.cpp
$(WRK)/test/benchmark/main.o: $(SRC)/test/benchmark/main.cpp; $(CXX_RULE) $(SRC)/test/benchmark/main.cpp
$(WRK)/test/benchmark_mpi/main.o: $(SRC)/test/benchmark_mpi/main.cpp; $(CXX_RULE) $(SRC)/test/benchmark_mpi/main.cpp
$(WRK)/test/benchmark_mpiraw/main.o: $(SRC)/test/benchmark_mpiraw/main.cpp; $(CXX_RULE) $(SRC)/test/benchmark_mpiraw/main.cpp
$(WRK)/test/benchmark_normal/main.o: $(SRC)/test/benchmark_normal/main.cpp; $(CXX_RULE) $(SRC)/test/benchmark_normal/main.cpp
$(WRK)/test/color/main.o: $(SRC)/test/color/main.cpp; $(CXX_RULE) $(SRC)/test/color/main.cpp
$(WRK)/test/commhalo/main.o: $(SRC)/test/commhalo/main.cpp; $(CXX_RULE) $(SRC)/test/commhalo/main.cpp
$(WRK)/test/comm/main.o: $(SRC)/test/comm/main.cpp; $(CXX_RULE) $(SRC)/test/comm/main.cpp
$(WRK)/test/commmap/main.o: $(SRC)/test/commmap/main.cpp; $(CXX_RULE) $(SRC)/test/commmap/main.cpp
$(WRK)/test/commmap/manager.o: $(SRC)/test/commmap/manager.cpp; $(CXX_RULE) $(SRC)/test/commmap/manager.cpp
$(WRK)/test/commmap/rank.o: $(SRC)/test/commmap/rank.cpp; $(CXX_RULE) $(SRC)/test/commmap/rank.cpp
$(WRK)/test/condface/main.o: $(SRC)/test/condface/main.cpp; $(CXX_RULE) $(SRC)/test/condface/main.cpp
$(WRK)/test/debug/main.o: $(SRC)/test/debug/main.cpp; $(CXX_RULE) $(SRC)/test/debug/main.cpp
$(WRK)/test/dump/dump_diff.o: $(SRC)/test/dump/dump_diff.cpp; $(CXX_RULE) $(SRC)/test/dump/dump_diff.cpp
$(WRK)/test/dump/dump_gen.o: $(SRC)/test/dump/dump_gen.cpp; $(CXX_RULE) $(SRC)/test/dump/dump_gen.cpp
$(WRK)/test/dump/dump_meta.o: $(SRC)/test/dump/dump_meta.cpp; $(CXX_RULE) $(SRC)/test/dump/dump_meta.cpp
$(WRK)/test/dump/dump_util.o: $(SRC)/test/dump/dump_util.cpp; $(CXX_RULE) $(SRC)/test/dump/dump_util.cpp
$(WRK)/test/embed_approx/main.o: $(SRC)/test/embed_approx/main.cpp; $(CXX_RULE) $(SRC)/test/embed_approx/main.cpp
$(WRK)/test/embed_convdiff/main.o: $(SRC)/test/embed_convdiff/main.cpp; $(CXX_RULE) $(SRC)/test/embed_convdiff/main.cpp
$(WRK)/test/embed_grad/main.o: $(SRC)/test/embed_grad/main.cpp; $(CXX_RULE) $(SRC)/test/embed_grad/main.cpp
$(WRK)/test/embed_interpolate/main.o: $(SRC)/test/embed_interpolate/main.cpp; $(CXX_RULE) $(SRC)/test/embed_interpolate/main.cpp
$(WRK)/test/embed/main.o: $(SRC)/test/embed/main.cpp; $(CXX_RULE) $(SRC)/test/embed/main.cpp
$(WRK)/test/embed_potential/main.o: $(SRC)/test/embed_potential/main.cpp; $(CXX_RULE) $(SRC)/test/embed_potential/main.cpp
$(WRK)/test/fluid/main.o: $(SRC)/test/fluid/main.cpp; $(CXX_RULE) $(SRC)/test/fluid/main.cpp
$(WRK)/test/grad/main.o: $(SRC)/test/grad/main.cpp; $(CXX_RULE) $(SRC)/test/grad/main.cpp
$(WRK)/test/hook/hook/hook.o: $(SRC)/test/hook/hook/hook.cpp; $(CXX_RULE) $(SRC)/test/hook/hook/hook.cpp
$(WRK)/test/hypre/main.o: $(SRC)/test/hypre/main.cpp; $(CXX_RULE) $(SRC)/test/hypre/main.cpp
$(WRK)/test/hypresub/main.o: $(SRC)/test/hypresub/main.cpp; $(CXX_RULE) $(SRC)/test/hypresub/main.cpp
$(WRK)/test/init/main.o: $(SRC)/test/init/main.cpp; $(CXX_RULE) $(SRC)/test/init/main.cpp
$(WRK)/test/inner/main.o: $(SRC)/test/inner/main.cpp; $(CXX_RULE) $(SRC)/test/inner/main.cpp
$(WRK)/test/linear/main.o: $(SRC)/test/linear/main.cpp; $(CXX_RULE) $(SRC)/test/linear/main.cpp
$(WRK)/test/mesh/main.o: $(SRC)/test/mesh/main.cpp; $(CXX_RULE) $(SRC)/test/mesh/main.cpp
$(WRK)/test/openmp/main.o: $(SRC)/test/openmp/main.cpp; $(CXX_RULE) $(SRC)/test/openmp/main.cpp
$(WRK)/test/parser/codeblocks.o: $(SRC)/test/parser/codeblocks.cpp; $(CXX_RULE) $(SRC)/test/parser/codeblocks.cpp
$(WRK)/test/parser/main.o: $(SRC)/test/parser/main.cpp; $(CXX_RULE) $(SRC)/test/parser/main.cpp
$(WRK)/test/particles/main.o: $(SRC)/test/particles/main.cpp; $(CXX_RULE) $(SRC)/test/particles/main.cpp
$(WRK)/test/primlist/getmap.o: $(SRC)/test/primlist/getmap.cpp; $(CXX_RULE) $(SRC)/test/primlist/getmap.cpp
$(WRK)/test/primlist/parse.o: $(SRC)/test/primlist/parse.cpp; $(CXX_RULE) $(SRC)/test/primlist/parse.cpp
$(WRK)/test/range/main.o: $(SRC)/test/range/main.cpp; $(CXX_RULE) $(SRC)/test/range/main.cpp
$(WRK)/test/recolor/main.o: $(SRC)/test/recolor/main.cpp; $(CXX_RULE) $(SRC)/test/recolor/main.cpp
$(WRK)/test/reconst/levelset.o: $(SRC)/test/reconst/levelset.cpp; $(CXX_RULE) $(SRC)/test/reconst/levelset.cpp
$(WRK)/test/reconst/main.o: $(SRC)/test/reconst/main.cpp; $(CXX_RULE) $(SRC)/test/reconst/main.cpp
$(WRK)/test/reconst/plane.o: $(SRC)/test/reconst/plane.cpp; $(CXX_RULE) $(SRC)/test/reconst/plane.cpp
$(WRK)/test/reflect/main.o: $(SRC)/test/reflect/main.cpp; $(CXX_RULE) $(SRC)/test/reflect/main.cpp
$(WRK)/test/stat/main.o: $(SRC)/test/stat/main.cpp; $(CXX_RULE) $(SRC)/test/stat/main.cpp
$(WRK)/test/subcomm/main.o: $(SRC)/test/subcomm/main.cpp; $(CXX_RULE) $(SRC)/test/subcomm/main.cpp
$(WRK)/test/suspender/main.o: $(SRC)/test/suspender/main.cpp; $(CXX_RULE) $(SRC)/test/suspender/main.cpp
$(WRK)/test/sysinfo/main.o: $(SRC)/test/sysinfo/main.cpp; $(CXX_RULE) $(SRC)/test/sysinfo/main.cpp
$(WRK)/test/sysinfo/realpath.o: $(SRC)/test/sysinfo/realpath.cpp; $(CXX_RULE) $(SRC)/test/sysinfo/realpath.cpp
$(WRK)/test/sysinfo/splitext.o: $(SRC)/test/sysinfo/splitext.cpp; $(CXX_RULE) $(SRC)/test/sysinfo/splitext.cpp
$(WRK)/test/tracker/main.o: $(SRC)/test/tracker/main.cpp; $(CXX_RULE) $(SRC)/test/tracker/main.cpp
$(WRK)/test/util/main.o: $(SRC)/test/util/main.cpp; $(CXX_RULE) $(SRC)/test/util/main.cpp
$(WRK)/test/vars/main.o: $(SRC)/test/vars/main.cpp; $(CXX_RULE) $(SRC)/test/vars/main.cpp
$(WRK)/test/visual/main.o: $(SRC)/test/visual/main.cpp; $(CXX_RULE) $(SRC)/test/visual/main.cpp
$(WRK)/util/convdiff.o: $(SRC)/util/convdiff.cpp; $(CXX_RULE) $(SRC)/util/convdiff.cpp
$(WRK)/util/distr.o: $(SRC)/util/distr.cpp; $(CXX_RULE) $(SRC)/util/distr.cpp
$(WRK)/util/events.o: $(SRC)/util/events.cpp; $(CXX_RULE) $(SRC)/util/events.cpp
$(WRK)/util/filesystem.o: $(SRC)/util/filesystem.cpp; $(CXX_RULE) $(SRC)/util/filesystem.cpp
$(WRK)/util/fixed_allocator.o: $(SRC)/util/fixed_allocator.cpp; $(CXX_RULE) $(SRC)/util/fixed_allocator.cpp
$(WRK)/util/fluid.o: $(SRC)/util/fluid.cpp; $(CXX_RULE) $(SRC)/util/fluid.cpp
$(WRK)/util/format.o: $(SRC)/util/format.cpp; $(CXX_RULE) $(SRC)/util/format.cpp
$(WRK)/util/gitgen.o: $(SRC)/util/gitgen.cpp; $(CXX_RULE) $(SRC)/util/gitgen.cpp
$(WRK)/util/git.o: $(SRC)/util/git.cpp; $(CXX_RULE) $(SRC)/util/git.cpp
$(WRK)/util/histogram.o: $(SRC)/util/histogram.cpp; $(CXX_RULE) $(SRC)/util/histogram.cpp
$(WRK)/util/hydro_post.o: $(SRC)/util/hydro_post.cpp; $(CXX_RULE) $(SRC)/util/hydro_post.cpp
$(WRK)/util/hydro.o: $(SRC)/util/hydro.cpp; $(CXX_RULE) $(SRC)/util/hydro.cpp
$(WRK)/util/linear.o: $(SRC)/util/linear.cpp; $(CXX_RULE) $(SRC)/util/linear.cpp
$(WRK)/util/logger.o: $(SRC)/util/logger.cpp; $(CXX_RULE) $(SRC)/util/logger.cpp
$(WRK)/util/mpi.o: $(SRC)/util/mpi.cpp; $(CXX_RULE) $(SRC)/util/mpi.cpp
$(WRK)/util/posthook_default.o: $(SRC)/util/posthook_default.cpp; $(CXX_RULE) $(SRC)/util/posthook_default.cpp
$(WRK)/util/suspender.o: $(SRC)/util/suspender.cpp; $(CXX_RULE) $(SRC)/util/suspender.cpp
$(WRK)/util/sysinfo.o: $(SRC)/util/sysinfo.cpp; $(CXX_RULE) $(SRC)/util/sysinfo.cpp
$(WRK)/util/system.o: $(SRC)/util/system.c; $(CC_RULE) $(SRC)/util/system.c
$(WRK)/util/timer.o: $(SRC)/util/timer.cpp; $(CXX_RULE) $(SRC)/util/timer.cpp
$(WRK)/util/visual.o: $(SRC)/util/visual.cpp; $(CXX_RULE) $(SRC)/util/visual.cpp
$(WRK)/util/vof.o: $(SRC)/util/vof.cpp; $(CXX_RULE) $(SRC)/util/vof.cpp
$(WRK)/young/main.o: $(SRC)/young/main.c; $(CC_RULE) $(SRC)/young/main.c
$(WRK)/explorer: $(WRK)/explorer.o; $(LINK) $(WRK)/explorer.o $(LINK_FLAGS)
$(WRK)/main: $(WRK)/main.o; $(LINK) $(WRK)/main.o $(LINK_FLAGS)
$(WRK)/test/advection/main: $(WRK)/test/advection/main.o; $(LINK) $(WRK)/test/advection/main.o $(LINK_FLAGS)
$(WRK)/test/approx/linear: $(WRK)/test/approx/linear.o; $(LINK) $(WRK)/test/approx/linear.o $(LINK_FLAGS)
$(WRK)/test/approx/main: $(WRK)/test/approx/main.o; $(LINK) $(WRK)/test/approx/main.o $(LINK_FLAGS)
$(WRK)/test/benchmark_diffusion/main: $(WRK)/test/benchmark_diffusion/main.o; $(LINK) $(WRK)/test/benchmark_diffusion/main.o $(LINK_FLAGS)
$(WRK)/test/benchmark/main: $(WRK)/test/benchmark/main.o; $(LINK) $(WRK)/test/benchmark/main.o $(LINK_FLAGS)
$(WRK)/test/benchmark_mpi/main: $(WRK)/test/benchmark_mpi/main.o; $(LINK) $(WRK)/test/benchmark_mpi/main.o $(LINK_FLAGS)
$(WRK)/test/benchmark_mpiraw/main: $(WRK)/test/benchmark_mpiraw/main.o; $(LINK) $(WRK)/test/benchmark_mpiraw/main.o $(LINK_FLAGS)
$(WRK)/test/benchmark_normal/main: $(WRK)/test/benchmark_normal/main.o; $(LINK) $(WRK)/test/benchmark_normal/main.o $(LINK_FLAGS)
$(WRK)/test/color/main: $(WRK)/test/color/main.o; $(LINK) $(WRK)/test/color/main.o $(LINK_FLAGS)
$(WRK)/test/commhalo/main: $(WRK)/test/commhalo/main.o; $(LINK) $(WRK)/test/commhalo/main.o $(LINK_FLAGS)
$(WRK)/test/comm/main: $(WRK)/test/comm/main.o; $(LINK) $(WRK)/test/comm/main.o $(LINK_FLAGS)
$(WRK)/test/commmap/main: $(WRK)/test/commmap/main.o; $(LINK) $(WRK)/test/commmap/main.o $(LINK_FLAGS)
$(WRK)/test/commmap/manager: $(WRK)/test/commmap/manager.o; $(LINK) $(WRK)/test/commmap/manager.o $(LINK_FLAGS)
$(WRK)/test/commmap/rank: $(WRK)/test/commmap/rank.o; $(LINK) $(WRK)/test/commmap/rank.o $(LINK_FLAGS)
$(WRK)/test/condface/main: $(WRK)/test/condface/main.o; $(LINK) $(WRK)/test/condface/main.o $(LINK_FLAGS)
$(WRK)/test/debug/main: $(WRK)/test/debug/main.o; $(LINK) $(WRK)/test/debug/main.o $(LINK_FLAGS)
$(WRK)/test/dump/dump_diff: $(WRK)/test/dump/dump_diff.o; $(LINK) $(WRK)/test/dump/dump_diff.o $(LINK_FLAGS)
$(WRK)/test/dump/dump_gen: $(WRK)/test/dump/dump_gen.o; $(LINK) $(WRK)/test/dump/dump_gen.o $(LINK_FLAGS)
$(WRK)/test/dump/dump_meta: $(WRK)/test/dump/dump_meta.o; $(LINK) $(WRK)/test/dump/dump_meta.o $(LINK_FLAGS)
$(WRK)/test/dump/dump_util: $(WRK)/test/dump/dump_util.o; $(LINK) $(WRK)/test/dump/dump_util.o $(LINK_FLAGS)
$(WRK)/test/embed_approx/main: $(WRK)/test/embed_approx/main.o; $(LINK) $(WRK)/test/embed_approx/main.o $(LINK_FLAGS)
$(WRK)/test/embed_convdiff/main: $(WRK)/test/embed_convdiff/main.o; $(LINK) $(WRK)/test/embed_convdiff/main.o $(LINK_FLAGS)
$(WRK)/test/embed_grad/main: $(WRK)/test/embed_grad/main.o; $(LINK) $(WRK)/test/embed_grad/main.o $(LINK_FLAGS)
$(WRK)/test/embed_interpolate/main: $(WRK)/test/embed_interpolate/main.o; $(LINK) $(WRK)/test/embed_interpolate/main.o $(LINK_FLAGS)
$(WRK)/test/embed/main: $(WRK)/test/embed/main.o; $(LINK) $(WRK)/test/embed/main.o $(LINK_FLAGS)
$(WRK)/test/embed_potential/main: $(WRK)/test/embed_potential/main.o; $(LINK) $(WRK)/test/embed_potential/main.o $(LINK_FLAGS)
$(WRK)/test/fluid/main: $(WRK)/test/fluid/main.o; $(LINK) $(WRK)/test/fluid/main.o $(LINK_FLAGS)
$(WRK)/test/grad/main: $(WRK)/test/grad/main.o; $(LINK) $(WRK)/test/grad/main.o $(LINK_FLAGS)
$(WRK)/test/hook/hook/hook: $(WRK)/test/hook/hook/hook.o; $(LINK) $(WRK)/test/hook/hook/hook.o $(LINK_FLAGS)
$(WRK)/test/hypre/main: $(WRK)/test/hypre/main.o; $(LINK) $(WRK)/test/hypre/main.o $(LINK_FLAGS)
$(WRK)/test/hypresub/main: $(WRK)/test/hypresub/main.o; $(LINK) $(WRK)/test/hypresub/main.o $(LINK_FLAGS)
$(WRK)/test/init/main: $(WRK)/test/init/main.o; $(LINK) $(WRK)/test/init/main.o $(LINK_FLAGS)
$(WRK)/test/inner/main: $(WRK)/test/inner/main.o; $(LINK) $(WRK)/test/inner/main.o $(LINK_FLAGS)
$(WRK)/test/linear/main: $(WRK)/test/linear/main.o; $(LINK) $(WRK)/test/linear/main.o $(LINK_FLAGS)
$(WRK)/test/mesh/main: $(WRK)/test/mesh/main.o; $(LINK) $(WRK)/test/mesh/main.o $(LINK_FLAGS)
$(WRK)/test/openmp/main: $(WRK)/test/openmp/main.o; $(LINK) $(WRK)/test/openmp/main.o $(LINK_FLAGS)
$(WRK)/test/parser/codeblocks: $(WRK)/test/parser/codeblocks.o; $(LINK) $(WRK)/test/parser/codeblocks.o $(LINK_FLAGS)
$(WRK)/test/parser/main: $(WRK)/test/parser/main.o; $(LINK) $(WRK)/test/parser/main.o $(LINK_FLAGS)
$(WRK)/test/particles/main: $(WRK)/test/particles/main.o; $(LINK) $(WRK)/test/particles/main.o $(LINK_FLAGS)
$(WRK)/test/primlist/getmap: $(WRK)/test/primlist/getmap.o; $(LINK) $(WRK)/test/primlist/getmap.o $(LINK_FLAGS)
$(WRK)/test/primlist/parse: $(WRK)/test/primlist/parse.o; $(LINK) $(WRK)/test/primlist/parse.o $(LINK_FLAGS)
$(WRK)/test/range/main: $(WRK)/test/range/main.o; $(LINK) $(WRK)/test/range/main.o $(LINK_FLAGS)
$(WRK)/test/recolor/main: $(WRK)/test/recolor/main.o; $(LINK) $(WRK)/test/recolor/main.o $(LINK_FLAGS)
$(WRK)/test/reconst/levelset: $(WRK)/test/reconst/levelset.o; $(LINK) $(WRK)/test/reconst/levelset.o $(LINK_FLAGS)
$(WRK)/test/reconst/main: $(WRK)/test/reconst/main.o; $(LINK) $(WRK)/test/reconst/main.o $(LINK_FLAGS)
$(WRK)/test/reconst/plane: $(WRK)/test/reconst/plane.o; $(LINK) $(WRK)/test/reconst/plane.o $(LINK_FLAGS)
$(WRK)/test/reflect/main: $(WRK)/test/reflect/main.o; $(LINK) $(WRK)/test/reflect/main.o $(LINK_FLAGS)
$(WRK)/test/stat/main: $(WRK)/test/stat/main.o; $(LINK) $(WRK)/test/stat/main.o $(LINK_FLAGS)
$(WRK)/test/subcomm/main: $(WRK)/test/subcomm/main.o; $(LINK) $(WRK)/test/subcomm/main.o $(LINK_FLAGS)
$(WRK)/test/suspender/main: $(WRK)/test/suspender/main.o; $(LINK) $(WRK)/test/suspender/main.o $(LINK_FLAGS)
$(WRK)/test/sysinfo/main: $(WRK)/test/sysinfo/main.o; $(LINK) $(WRK)/test/sysinfo/main.o $(LINK_FLAGS)
$(WRK)/test/sysinfo/realpath: $(WRK)/test/sysinfo/realpath.o; $(LINK) $(WRK)/test/sysinfo/realpath.o $(LINK_FLAGS)
$(WRK)/test/sysinfo/splitext: $(WRK)/test/sysinfo/splitext.o; $(LINK) $(WRK)/test/sysinfo/splitext.o $(LINK_FLAGS)
$(WRK)/test/tracker/main: $(WRK)/test/tracker/main.o; $(LINK) $(WRK)/test/tracker/main.o $(LINK_FLAGS)
$(WRK)/test/util/main: $(WRK)/test/util/main.o; $(LINK) $(WRK)/test/util/main.o $(LINK_FLAGS)
$(WRK)/test/vars/main: $(WRK)/test/vars/main.o; $(LINK) $(WRK)/test/vars/main.o $(LINK_FLAGS)
$(WRK)/test/visual/main: $(WRK)/test/visual/main.o; $(LINK) $(WRK)/test/visual/main.o $(LINK_FLAGS)
