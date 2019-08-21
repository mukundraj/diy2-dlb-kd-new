#!/bin/bash

#SBATCH --job-name=baseline
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=64
#SBATCH --output=baseline.%j.%N.out
#SBATCH --error=baseline.%j.%N.error
#SBATCH --mail-user=mraj@lcrc.anl.gov
#SBATCH --time=01:00:00

srun -n 2 ./sync_nek configs_48/config.nek2.xml 2> 48/nek_02
srun -n 4 ./sync_nek configs_48/config.nek4.xml 2> 48/nek_04
srun -n 8 ./sync_nek configs_48/config.nek8.xml 2> 48/nek_08
srun -n 16 ./sync_nek configs_48/config.nek16.xml 2> 48/nek_16
srun -n 32 ./sync_nek configs_48/config.nek32.xml 2> 48/nek_32
srun -n 64 ./sync_nek configs_48/config.nek64.xml 2> 48/nek_64
srun -n 128 ./sync_nek configs_48/config.nek128.xml 2> 48/nek_128
srun -n 256 ./sync_nek configs_48/config.nek256.xml 2> 48/nek_256
srun -n 512 ./sync_nek configs_48/config.nek512.xml 2> 48/nek_512

# srun -n 2 ./sync_nek configs_96/config.nek2.xml 2> 96/nek_02
# srun -n 4 ./sync_nek configs_96/config.nek4.xml 2> 96/nek_04
# srun -n 8 ./sync_nek configs_96/config.nek8.xml 2> 96/nek_08
# srun -n 16 ./sync_nek configs_96/config.nek16.xml 2> 96/nek_16
# srun -n 32 ./sync_nek configs_96/config.nek32.xml 2> 96/nek_32
# srun -n 64 ./sync_nek configs_96/config.nek64.xml 2> 96/nek_64
# srun -n 128 ./sync_nek configs_96/config.nek128.xml 2> 96/nek_128
# srun -n 256 ./sync_nek configs_96/config.nek256.xml 2> 96/nek_256
# srun -n 512 ./sync_nek configs_96/config.nek512.xml 2> 96/nek_512

# srun -n 2 ./sync_nek configs_384/config.nek2.xml 2> 128/nek_02
# srun -n 4 ./sync_nek configs_384/config.nek4.xml 2> 128/nek_04
# srun -n 8 ./sync_nek configs_384/config.nek8.xml 2> 128/nek_08
# srun -n 16 ./sync_nek configs_384/config.nek16.xml 2> 128/nek_16
# srun -n 32 ./sync_nek configs_384/config.nek32.xml 2> 128/nek_32
# srun -n 64 ./sync_nek configs_384/config.nek64.xml 2> 128/nek_64
# srun -n 128 ./sync_nek configs_384/config.nek128.xml 2> 128/nek_128
# srun -n 256 ./sync_nek configs_384/config.nek256.xml 2> 128/nek_256
# srun -n 512 ./sync_nek configs_384/config.nek512.xml 2> 128/nek_512