#import matplotlib.pyplot as p
from sys import argv
from textwrap import fill
from glob import glob

domainmin=40
scorethresh=1.4
lenlimit=1500
outdir=argv[1]

allprots=glob("%s/fasta/*/"%outdir)
for prot in allprots:
    protname=prot.split("/")[-2]
    g=None
    try:
        g=open("%s/seq.ConDo"%prot)
    except:
        #print("ConDo did not run properly for target %s"%prot)
        fasta=open("%s/seq.fasta"%prot)
        fasta.readline()
        seq=fasta.readline().strip()
        fasta.close()
        if len(seq) > lenlimit:
            print("Warning: Target %s is longer than SPRING can handle.  Only the first %d amino acids were written."%(protname,lenlimit))
        outf=open("%s/1.fasta"%prot,"w")
        outf.write(">%s|1\n"%protname)
        outf.write(fill(seq[:lenlimit],width=60))
        outf.close()
        continue
    vals=[float(line.split()[1]) for line in g]
    g.close()
    boundaries=[]
    highestind=-1
    lastind=-1
    highestval=scorethresh
    for i in range(len(vals)):
        if vals[i]>scorethresh:
            lastind=i
            if vals[i]>highestval:
                highestind=i
                highestval=vals[i]
        else:
            if highestind>=0 and i-lastind>=domainmin:
                boundaries.append(highestind)
                highestind=-1
                lastind=-1
                highestval=scorethresh
    fasta=open("%s/seq.fasta"%prot)
    fasta.readline()
    seq=fasta.readline().strip()
    fasta.close()
    boundaries=[0]+boundaries+[len(seq)]
    for i in range(len(boundaries)-1):
        outf=open("%s/%d.fasta"%(prot,i+1),"w")
        outf.write(">%s|%d\n"%(protname,i+1))
        seqslice=seq[boundaries[i]:boundaries[i+1]]
        if len(seqslice) > lenlimit:
            print("Warning: Target %s_%d is longer than SPRING can handle.  Only the first %d amino acids were written."%(protname,i+1,lenlimit))
        outf.write(fill(seqslice[:lenlimit],width=60))
        outf.close()
