#!/bin/bash

outdir=$1
prog=$2

cat ${outdir}/PPI/${prog}res_*.txt >> ${outdir}/PPI/${prog}res.txt
rm -rf ${outdir}/PPI/${prog}res_*.txt
sort ${outdir}/PPI/${prog}res.txt | uniq > ${outdir}/PPI/m
mv ${outdir}/PPI/m ${outdir}/PPI/${prog}res.txt
