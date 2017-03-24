#!/bin/bash

# Chaste template script
# Use this template to create a script for running tests on supercomputers
#  * First edit this file to suit your system and job scheduler
#  * Then compile on the head node with this command: "scons offline_mode=1 <list of tests>"

# EDIT HERE - commands for the scheduler
# Here is an example to run in the current directory (works on GridEngine)
#$ -cwd
# GridEngine join output and error files
#$ -j oe
# PBS join output and error files
#PBS -j oe

# EDIT HERE - launch command
# Define the command and any other machine files
# Here is an example for some supercomputer from 7 years ago:
# export MPI_LAUNCH_COMMAND="mpirun -np $NSLOTS -machinefile $TMPDIR/machines"
export MPI_LAUNCH_COMMAND="mpirun"

# EDIT HERE - any other miscellaneous commands that need to be run before the tests
# For example in PBS:
# cd $PBS_O_WORKDIR


export LD_LIBRARY_PATH=/home/u1221/release_3.4/lib:/opt/chaste/petsc-3.4.4/linux-gnu/lib:/opt/chaste/lib:/usr/lib64/vtk:/home/u1221/release_3.4/lib
export PATH=.:/opt/chaste/petsc-3.4.4/linux-gnu/bin:/opt/intel/compilers_and_libraries_2016.0.109/linux/bin/intel64:/opt/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/bin:/opt/intel/debugger_2016/gdb/intel64_mic/bin:/opt/chaste/petsc-3.4.4/linux-gnu/bin:/usr/lib64/openmpi/bin:/opt/intel/compilers_and_libraries_2016.0.109/linux/bin/intel64:/opt/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/bin:/opt/intel/debugger_2016/gdb/intel64_mic/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/u1221/release_3.4/apps/src:/home/u1221/.local/bin:/home/u1221/bin:/home/u1221/release_3.4/apps/src


sbatch -t 2-23:00:00 --exclusive -N 1 --wrap "$MPI_LAUNCH_COMMAND heart/build/optimised_ndebug/nastya/nastya_lrRunner"
sleep 1
