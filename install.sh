#!/bin/bash

echo "Thank you for installing PEPPI! This installation will take ~100GB of space due to database size."
read -p "Are you on a slurm HPC system? (WARNING: PEPPI will run slowly without HPC parallelization) [y/n] " hpcanswer
if [ $hpcanswer = "n" ] || [ $hpcanswer = "no" ] || [ $hpcanswer = "N" ] || [ $hpcanswer = "No" ]; then
    hpcflag=0
    maxjobs=1
elif [ $hpcanswer = "y" ] || [ $hpcanswer = "yes" ] || [ $hpcanswer = "Y" ] || [ $hpcanswer = "Yes" ]; then
    hpcflag=1
    read -p "How many jobs are you able to run at once? " maxjobs 
else
    echo "Invalid response: $hpcanswer"
    exit 1
fi

echo ""
read -p "Full path of where you wish to install PEPPI: " peppidir
read -p "Full path to your HHsuite installation: " hhdir
read -p "Full path to the database used for hhblits: " dbdir
read -p "Full path of your python interpreter: " pythonbin
read -p "What is your C++ compiler? " cppcompiler
read -p "What is your fortran compiler? " fcompiler

mkdir -p $peppidir
cd $peppidir
if ! [ -d "PEPPI" ]; then
    git clone https://github.com/ewbell94/PEPPI.git
else
    echo "PEPPI source already exists!"
fi

cd PEPPI
if ! [ -d "lib" ]; then
    if ! [ -f "lib.tar.gz" ]; then
	wget https://zhanggroup.org/PEPPI/download/lib.tar.gz
    fi
    tar -zxvf lib.tar.gz && rm -rf lib.tar.gz
else
    echo "PEPPI library already exists!"
fi
sed -i "s#\$peppidir= \".*\"#\$peppidir=\"$peppidir/PEPPI\"#" PEPPI1.pl
sed -i "s#\$maxjobs=.*;#\$maxjobs=$maxjobs;#" PEPPI1.pl
sed -i "s#\$hpcflag=.*;#\$hpcflag=$hpcflag;#" PEPPI1.pl
sed -i "s#\$hhsuitedir=\".*\"#\$hhsuitedir=\"$hhdir\"#" bin/makeHHR.pl
sed -i "s#\$uniprotdb=\".*\"#\$uniprotdb=\"$dbdir\"#" bin/makeHHR.pl
sed -i "s#/nfs/amino-library/anaconda/bin/python#$pythonbin#" bin/CTmod
sed -i "s#/nfs/amino-library/anaconda/bin/python#$pythonbin#" bin/STRINGmod
sed -i "s#/nfs/amino-library/anaconda/bin/python#$pythonbin#" bin/PEPPI3temp.py
sed -i "s#/nfs/amino-home/ewbell#$peppidir#g" bin/compiled_source/dcomplex.c
$cppcompiler bin/compiled_source/dcomplex.c -o bin/dcomplex -lm -O3
$cppcompiler bin/compiled_source/dimMap.cpp -o bin/dimMap -O3 --std=c++11
$fcompiler bin/compiled_source/NWalign.f -o bin/NWalign -lm -O3
$fcompiler bin/compiled_source/TMalign.f -o bin/TMalign -lm -O3

rm -rf bin/model_multiD
rm -rf bin/CTNN
$pythonbin bin/trainCT.py
$pythonbin bin/trainDists.py

test=`$pythonbin bin/getHashcode.py "Testing" "Hashcode"`
if [ $test != "636" ]; then
    echo "WARNING: getHashcode.py is not functioning properly, this will cause STRING to malfunction."
fi
