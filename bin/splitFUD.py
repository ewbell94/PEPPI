from sys import argv
from textwrap import fill
from glob import glob

lenlimit=10000
outdir=argv[1]
allprots=glob("%s/fasta/*/"%outdir)

for prot in allprots:
    target=prot.split("/")[-2]    
    try:
        f=open(prot+"fu.txt")
        fures=f.readline()
        f.close()
    except:
        print("Target %s has no FUpred result"%target)
        fasta=open(prot+"seq.fasta")
        fasta.readline()
        seq=""
        for line in fasta:
            seq+=line.strip()
        if len(seq) > lenlimit:
            print("Warning: Target %s is longer than PEPPI can handle.  Only the first %d amino acids were written."%(target,lenlimit))
        outf=open(prot+"1.fasta","w")
        outf.write(">%s|1\n"%target)
        outf.write(fill(seq[:lenlimit]))
        outf.close()
        continue

    domains=fures.split(";")[:-1]
    f=open(prot+"seq.fasta")
    f.readline()
    seq=""
    for line in f:
        seq+=line.strip()
    for i in range(len(domains)):
        g=open(prot+"%d.fasta"%i,"w")
        g.write(">%s|%d\n"%(target,i))
        domainseq=""
        for part in domains[i].split(","):
            bounds=[int(j) for j in part.split("-")]
            domainseq+=seq[bounds[0]-1:bounds[1]]
        if len(domainseq) > lenlimit:
            print("Warning: Domain %s_%d is longer than PEPPI can handle.  Only the first %d amino acids were written."%(target,i,lenlimit))
        g.write(fill(domainseq[:lenlimit]))
        g.close()



