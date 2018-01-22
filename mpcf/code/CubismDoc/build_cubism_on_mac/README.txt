5 steps to run CUBISM-MPCF on your MAC:
===========================================

1. checkout out a clean version of CUBISM-MPCF from the git repository

2. apply the included patch file: 
   patch -p1 < necessary_diffs_for_mac.patch

3. rebuild 'libfpzip.a':
   cd CubismApps/tools/fpzip/src
   make clean
   make

4. build 'mpcf-cluster':
   cd ../../MPCFcluster/makefile
   make cleanall
   make
   ALTERNATIVE OPTION (which however seems not to be supported by all tests)
   build 'mpcf-node':
   cd ../../MPCFnode/makefile
   make cleanall
   make
   In this case, you do not include MPI.

Note: these changes have to be reverted before committing into the git repository by using
=====

./prepare_commit.sh

(If you copy this file to another directory, you have to adapt the path therein.)

5. actually run CUBISM-MPCF
   ./run_cubims_mpcf_on_mac.sh
  If you have built 'mpcf-node', replace 'mpcf-cluster'. Then you do not have to explictely set -xpesize 1 -ypesize 1 -zpesize 1 and may omit them.
