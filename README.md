This repository holds the code for the Pipeline for the Extraction of Predicted Protein-protein Interactions (PEPPI).

# Prerequisites
Note: version numbers are provided to ensure proper functionality.  Other version configurations may also be compatible, but the provided ones are guaranteed to work replicably.

1. A python 2.7.16 environment (preferred Anaconda) with the following packages installed:
   + numpy 1.14.1
   + scipy 1.0.0
   + scikit-learn 0.18.1
2. A perl 5 interpreter (v5.16.3)
3. A C/C++ compiler compatible with the C++11 standard (g++ 4.8.5)
4. A fortran compiler (gfortran 4.8.5)
5. An installation of HHsuite3 (you will need to install psipred, blast, and a database as part of this, otherwise HHsearch will not work properly)

# Installation
Run the included install.sh script.  You will be asked several questions regarding the configuration of PEPPI, including install locations.  If you are not running PEPPI on an HPC environment, this script will create aliases for the "sbatch" and "squeue" commands, which will allow PEPPI to run properly.  After installing, please re-open your terminal and run these aliases before proceeding to ensure they're properly defined (squeue is echo, sbatch is the perl interpreter).