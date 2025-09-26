#!/bin/bash -x
#SBATCH --nodes=4
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=20
#SBATCH --output=./slurm_dales-%j.out
#SBATCH --error=./slurm_error-%j.out
#SBATCH --time=24:00:00
#SBATCH --account=hr-afc
#SBATCH --partition=batch


# number of nodes in $SLURM_NNODES (default: 1)
# number of tasks in $SLURM_NTASKS (default: 1)
# number of tasks per node in $SLURM_NTASKS_PER_NODE (default: 1)
# number of threads per task in $SLURM_CPUS_PER_TASK (default: 1)

echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS



# load modules
module load GCC  ParaStationMPI
module load netCDF-Fortran
module load HDF
module load netCDF
module load CMake 


srun --cpu-bind=map_cpu:0,1,2,3,4,5,6,7,8,9,24,25,26,27,28,29,30,31,32,33 ./dales4 > stdlog.txt
python3 change_nml.py

