# MISA-AKMC
A GPU implementation of atomic kinetic Monte Carlo (AKMC) method.
# Hardware
All of the experiments were conducted on a computer cluster with more than 4,000 compute nodes, each containing a 32-core AMD zen processor, four AMD Instinct M160 GPUs, and 128 GB of host memory. 

The processor has four Core Complex Dies, each configured as a NUMA node. Each GPU has 64 CUs and 16 GB of global memory. 

The high-speed computing network utilizes the HDR InfiniBand network, with a network speed of 200Gbps between computing nodes.

The operating system is NFS Server 3.2. The slurm version is 20.11.8-1.4.1-20220823. The version of Module greater than 3.2.10 is recommended.

# Software
MISA-AKMC relies on a number of open source libraries such as kiwi, googletest, xoshiro, libcomm, etc. 
You can use our pkg dependency management tool to download the dependency package or import the corresponding dependency package directly into the MISA-AKMC source vendor directory. 

The URL of the pkg package management tool software is: https://github.com/genshen/pkg

Generally, you can select the latest version.

#Datasets / Inputs
MISA-AKMC use yaml format file as configuration file. 
The program reads the data in the config.yaml file to initialize the environment that needs to be computed. 
Users can modify the relevant parameters in the config.yaml file to achieve the situation they want to simulate. 
A sample config.yaml file is available in the examples folder of the MISA-AKMC package.

The URL of an example of the config.yaml file in MISA-AKMC: https://github.com/panzhija/MISA-AKMC-GPU/blob/main/examples/config.yaml

# Installation and Deployment
Some of the necessary compilers we usually use are:

devtoolset: 7.3.1

intel: 2017.5.239

intelmpi: 2017.4.239

rocm/dtk: 23.04

cmake: 3.24.1

You can also use other compilers, but make sure they meet the following conditions:

CMake: 3.15 and above.

C++ compiler that supports C++11 features.

MPI Environment: An MPI environment that supports MPI 2.0 and later is required to be installed on your system.

Before the program can be compiled, we need to use the pkg C/C++ dependency management tool to process the required dependency packages. 
You can download it at the URL we provided above. For Linux operating systems with 64-bit amd64 architecture, we can ensure that the pkg executable is in the environment variable by entering the following command.

mkdir -p ~/.local/bin

wget https://github.com/genshen/pkg/releases/download/v0.4.1/pkg-linux-amd64 -O ~/.local/bin/pkg

chmod +x ~/.local/bin/pkg

export PATH=~/.local/bin:\$PATH

We can download MISA-AKMC from github or use "git clone https://github.com/panzhija/MISA-AKMC-GPU.git" command. 
Then, after "cd misa-akmc", we can use the "pkg fetch" and "pkg install" commands in turn to import the related dependency packages.
Finally, we can use "cmake -B./ build-H./" and "cmake --build./build" commands, the compilation is complete. 
MISA-AKMC is calculated by default using the GPU, if you want to use the CPU, you can use the command "cmake -B ./build -H./ -DKMC\_HIP\_ENABLE\_FLAG=OFF".

Experimental tasks are submitted via slurm scripts. An example of a slurm script is shown below.

\$ cat slurm.job

\#!/bin/bash

\#SBATCH -J misa-akmc

\#SBATCH -p noraml

\#SBATCH -N 2

\#SBATCH --exclusive

\#SBATCH --ntasks-per-node=4

\#SBATCH --mem=100G

\#SBATCH -o out.\%j

\#SBATCH -e \%j.loop

\#SBATCH --gres=dcu:4

echo "SLURM\_JOB\_NODELIST=\${SLURM\_JOB\_NODELIST}"

echo "SLURM\_NODELIST=\${SLURM\_NODELIST}"

EXEC\_PATH=/misa-akmc/build

CONFIG\_PATH=/test\_dcu/misa-akmc

mpirun \$EXEC\_PATH/bin/pmisa-kmc -c \$CONFIG\_PATH/config.yaml

date

The job can be submitted using the "sbatch slurm.job" command.

# Outputs
OVITO is a scientific visualization and analysis software for atomic and particle simulation data. 

It can be downloaded at "https://www.ovito.org/\#download".

The kmc.out file output by the MISA-AKMC contains the coordinates and atomic type of each lattice point, and one row of data is (atom\_type, x, y, z). 
We load the kmc.out file into the OVITO software and can see the visualized image of the simulation. 
In addition, other output information of MISA-AKMC can also be viewed in the result file, such as: the number of discrete solute atoms, program execution time, etc.
