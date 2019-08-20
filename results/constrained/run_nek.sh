#!/bin/sh

mpirun -np 2 ./sync_nek configs_512/config.nek2.xml 2> 512/nek_02
mpirun -np 4 ./sync_nek configs_512/config.nek4.xml 2> 512/nek_04
mpirun -np 8 ./sync_nek configs_512/config.nek8.xml 2> 512/nek_08
mpirun -np 16 ./sync_nek configs_512/config.nek16.xml 2> 512/nek_16
mpirun -np 32 ./sync_nek configs_512/config.nek32.xml 2> 512/nek_32

echo "half done"

mpirun -np 2 ./sync_nek configs_1024/config.nek2.xml 2> 1024/nek_02
mpirun -np 4 ./sync_nek configs_1024/config.nek4.xml 2> 1024/nek_04
mpirun -np 8 ./sync_nek configs_1024/config.nek8.xml 2> 1024/nek_08
mpirun -np 16 ./sync_nek configs_1024/config.nek16.xml 2> 1024/nek_16
mpirun -np 32 ./sync_nek configs_1024/config.nek32.xml 2> 1024/nek_32