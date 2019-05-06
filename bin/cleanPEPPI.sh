#!/bin/bash

peppidir="/nfs/amino-home/ewbell/PEPPI"
cd $peppidir

rm -rf fasta hhr SPRING PEPPI2.pl PEPPI3.pl protcode.csv res.csv
find . -name "*~" -delete