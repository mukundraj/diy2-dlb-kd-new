#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --output=p1.%j.%N.out
#SBATCH --error=p1.%j.%N.error
#SBATCH --mail-user=mraj@lcrc.anl.gov
#SBATCH --time=03:00:00


# srun -n 1 ./sync_nek configs_48/config.nek2.xml 2> 48/nek_01
# srun -n 2 ./sync_nek configs_48/config.nek2.xml 2> 48/nek_02
# srun -n 4 ./sync_nek configs_48/config.nek4.xml 2> 48/nek_04
# srun -n 8 ./sync_nek configs_48/config.nek8.xml 2> 48/nek_08
# srun -n 16 ./sync_nek configs_48/config.nek16.xml 2> 48/nek_16
# srun -n 32 ./sync_nek configs_48/config.nek32.xml 2> 48/nek_32
srun -n 64 ./sync_nek configs_48/config.nek64.xml 2> 48/nek_64
srun -n 128 ./sync_nek configs_48/config.nek128.xml 2> 48/nek_128
srun -n 256 ./sync_nek configs_48/config.nek256.xml 2> 48/nek_256
srun -n 512 ./sync_nek configs_48/config.nek512.xml 2> 48/nek_512
# srun -n 1024 ./sync_nek configs_48/config.nek1024.xml 2> 48/nek_1024
# srun -n 2048 ./sync_nek configs_48/config.nek2048.xml 2> 48/nek_2048
# srun -n 4096 ./sync_nek configs_48/config.nek4096.xml 2> 48/nek_4096

# srun -n 1 ./sync_nek configs_96/config.nek2.xml 2> 96/nek_01
# srun -n 2 ./sync_nek configs_96/config.nek2.xml 2> 96/nek_02
# srun -n 4 ./sync_nek configs_96/config.nek4.xml 2> 96/nek_04
# srun -n 8 ./sync_nek configs_96/config.nek8.xml 2> 96/nek_08
# srun -n 16 ./sync_nek configs_96/config.nek16.xml 2> 96/nek_16
# srun -n 32 ./sync_nek configs_96/config.nek32.xml 2> 96/nek_32
srun -n 64 ./sync_nek configs_96/config.nek64.xml 2> 96/nek_64
srun -n 128 ./sync_nek configs_96/config.nek128.xml 2> 96/nek_128
srun -n 256 ./sync_nek configs_96/config.nek256.xml 2> 96/nek_256
srun -n 512 ./sync_nek configs_96/config.nek512.xml 2> 96/nek_512
# srun -n 1024 ./sync_nek configs_96/config.nek1024.xml 2> 96/nek_1024
# srun -n 2048 ./sync_nek configs_96/config.nek2048.xml 2> 96/nek_2048
# srun -n 4096 ./sync_nek configs_96/config.nek4096.xml 2> 96/nek_4096


# srun -n 1 ./sync_nek configs_384/config.nek2.xml 2> 384/nek_01
# srun -n 2 ./sync_nek configs_384/config.nek2.xml 2> 384/nek_02
# srun -n 4 ./sync_nek configs_384/config.nek4.xml 2> 384/nek_04
# srun -n 8 ./sync_nek configs_384/config.nek8.xml 2> 384/nek_08
# srun -n 16 ./sync_nek configs_384/config.nek16.xml 2> 384/nek_16
# srun -n 32 ./sync_nek configs_384/config.nek32.xml 2> 384/nek_32
srun -n 64 ./sync_nek configs_384/config.nek64.xml 2> 384/nek_64
srun -n 128 ./sync_nek configs_384/config.nek128.xml 2> 384/nek_128
srun -n 256 ./sync_nek configs_384/config.nek256.xml 2> 384/nek_256
srun -n 512 ./sync_nek configs_384/config.nek512.xml 2> 384/nek_512
# srun -n 1024 ./sync_nek configs_384/config.nek1024.xml 2> 384/nek_1024
# srun -n 2048 ./sync_nek configs_384/config.nek2048.xml 2> 384/nek_2048
# srun -n 4096 ./sync_nek configs_384/config.nek4096.xml 2> 384/nek_4096


# srun -n 1 ./sync_nek configs_unlim/config.nek2.xml 2> unlim/nek_01
# srun -n 2 ./sync_nek configs_unlim/config.nek2.xml 2> unlim/nek_02
# srun -n 4 ./sync_nek configs_unlim/config.nek4.xml 2> unlim/nek_04
# srun -n 8 ./sync_nek configs_unlim/config.nek8.xml 2> unlim/nek_08
# srun -n 16 ./sync_nek configs_unlim/config.nek16.xml 2> unlim/nek_16
# srun -n 32 ./sync_nek configs_unlim/config.nek32.xml 2> unlim/nek_32
srun -n 64 ./sync_nek configs_unlim/config.nek64.xml 2> unlim/nek_64
srun -n 128 ./sync_nek configs_unlim/config.nek128.xml 2> unlim/nek_128
srun -n 256 ./sync_nek configs_unlim/config.nek256.xml 2> unlim/nek_256
srun -n 512 ./sync_nek configs_unlim/config.nek512.xml 2> unlim/nek_512
# srun -n 1024 ./sync_nek configs_unlim/config.nek1024.xml 2> unlim/nek_1024
# srun -n 2048 ./sync_nek configs_unlim/config.nek2048.xml 2> unlim/nek_2048
# srun -n 4096 ./sync_nek configs_unlim/config.nek4096.xml 2> unlim/nek_4096