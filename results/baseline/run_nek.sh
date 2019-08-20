#!/bin/sh

mpirun -n 2 ./sync_nek configs_512/config.nek2.xml 2> 512/nek_02
mpirun -n 4 ./sync_nek configs_512/config.nek4.xml 2> 512/nek_04
mpirun -n 8 ./sync_nek configs_512/config.nek8.xml 2> 512/nek_08
mpirun -n 16 ./sync_nek configs_512/config.nek16.xml 2> 512/nek_16
mpirun -n 32 ./sync_nek configs_512/config.nek32.xml 2> 512/nek_32
