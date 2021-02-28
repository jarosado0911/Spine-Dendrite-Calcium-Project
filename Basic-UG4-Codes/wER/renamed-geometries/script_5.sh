#!/bin/sh 
#PBS -l walltime=00:60:00 
#PBS -N Cell_5_RyR
#PBS -q normal 
#PBS -l nodes=1:ppn=28 
#PBS -e error.txt 
#PBS -o output.txt 
cd $PBS_O_WORKDIR 

mpirun -np 28 ugshell -ex /home/tuh36651/work/ActiveER-FixedCA/reconstructed_spine_wER.lua -grid /home/tuh36651/work/ActiveER-FixedCA/Spine5_wER.ugx -numRefs 0 -caInflux 0.0086 -tstep 5e-6 -endTime 0.030 -outName /home/tuh36651/work/ActiveER-FixedCA/output5 -solver GS -setting ryr -ryrDensity 3.13 -minDef 1e-12 -freq 1