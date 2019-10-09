#!/bin/bash -l
#SBATCH --job-name="wave_vofm2"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint=mc
#SBATCH --hint=multithread
#SBATCH --cpus-per-task=72
#SBATCH --ntasks-per-core=2
#SBATCH --exclusive
#SBATCH --time=5:00:00
##SBATCH --output=./out.log
##SBATCH --error=./err.log

module load daint-mc
module unload ddt
module use /scratch/snx3000/jfavre/daint/modules/all
module load ParaView/5.7.0-CrayGNU-18.08-OSMesa


echo "PYTHONPATH=" $PYTHONPATH
echo "----------------------------------------------------------------"
echo "LD_LIBRARY_PATH=" $LD_LIBRARY_PATH
echo "----------------------------------------------------------------"
echo "PATH=" $PATH
echo "----------------------------------------------------------------"

export VISRTX_TONE_MAPPING=0
srun -n $SLURM_NTASKS -N $SLURM_NNODES  `which pvbatch`  ./rtblue.py -fine2 ../sm_*.vtk 

