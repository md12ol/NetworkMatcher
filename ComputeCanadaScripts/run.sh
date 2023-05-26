#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=1G
#SBATCH --account=def-houghten
./cmake-build-release---remote/NetMatch 250 100000 20 3 29 1140