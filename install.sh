#!/bin/bash

read -p "Where do you wish to install PEPPI? " peppidir
read -p "Where is your HHsuite installation? " hhdir
read -p "How many jobs are you able to run at once? " maxjobs
read -p "What is your C++ compiler? " cppcompiler
read -p "What is your fortran compiler? " fcompiler

mkdir -p $peppidir
cd $peppidir
git clone https://github.com/ewbell94/PEPPI.git
cd PEPPI
wget https://zhanglab.dcmb.med.umich.edu/PEPPI/lib.tar.gz
tar -zxvf lib.tar.gz
sed -i "s#\$peppidir= \".*\"#\$peppidir=\"$peppidir\"#" PEPPI1.pl
sed -i "s#\$maxjobs=.*;#\$maxjobs=$maxjobs;#" PEPPI1.pl
sed -i "s#\$hhsuitedir=\".*\"#\$hhsuitedir=$hhdir#" bin/makeHHR.pl
$cppcompiler bin/compiled_source/dcomplex.c -o bin/dcomplex -lm -O3
$cppcompiler bin/compiled_source/dimMap.cpp -o bin/dimMap -O3 --std=c++11
$fcompiler bin/compiled_source/NWalign.f -o bin/NWalign -lm -O3
$fcompiler bin/compiled_source/TMalign.f -o bin/TMalign -lm -O3
