#!/bin/bash

echo "Thank you for installing PEPPI! This installation will take ~100GB of space due to database size."
read -p "Are you on a slurm HPC system? (WARNING: if no, this will overwrite sbatch and squeue) [y/n] " hpcflag
if [ $hpcflag = "n" ] || [ $hpcflag = "no" ]; then
    if [ $SHELL = "/bin/zsh" ]; then
	profile=~/.zshenv
    elif [ $SHELL = "/bin/csh" ] || [ $SHELL = "/bin/tcsh" ]; then
	profile=~/.login
    elif [ $SHELL = "/bin/bash" ]; then
	profile=~/.bashrc
    else
	profile=~/.profile
    fi
    echo "#The following changes were made by PEPPI" >> $profile
    echo "alias squeue='echo'" >> $profile
    echo "sbatch () {
    while [ \`echo \$1 | egrep -o \".{1,3}$\"\` != \".pl\" ]; do
        shift
    done
    perl \$@
}" >> $profile
    maxjobs=1

elif [ $hpcflag = "y" ] || [ $hpcflag = "yes" ]; then
    read -p "How many jobs are you able to run at once? " maxjobs 
else
    echo "Invalid response: $hpcflag"
    exit 1
fi

read -p "Where do you wish to install PEPPI? " peppidir
read -p "Where is your HHsuite installation? " hhdir
read -p "Where is the database used for hhblits? " dbdir
read -p "What is your C++ compiler? " cppcompiler
read -p "What is your fortran compiler? " fcompiler

mkdir -p $peppidir
cd $peppidir
wget https://zhanggroup.org/PEPPI/PEPPI.tar.gz
tar -zxvf PEPPI.tar.gz
cd PEPPI
wget https://zhanggroup.org/PEPPI/lib.tar.gz
tar -zxvf lib.tar.gz
sed -i "s#\$peppidir= \".*\"#\$peppidir=\"$peppidir\"#" PEPPI1.pl
sed -i "s#\$maxjobs=.*;#\$maxjobs=$maxjobs;#" PEPPI1.pl
sed -i "s#\$hhsuitedir=\".*\"#\$hhsuitedir=\"$hhdir\"#" bin/makeHHR.pl
sed -i "s#\$uniprotdb=\".*\"#\$uniprotdb=\"$dbdir\"" bin/makeHHR.pl
$cppcompiler bin/compiled_source/dcomplex.c -o bin/dcomplex -lm -O3
$cppcompiler bin/compiled_source/dimMap.cpp -o bin/dimMap -O3 --std=c++11
$fcompiler bin/compiled_source/NWalign.f -o bin/NWalign -lm -O3
$fcompiler bin/compiled_source/TMalign.f -o bin/TMalign -lm -O3

