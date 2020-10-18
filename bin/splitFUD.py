from sys import argv
from textwrap import fill
from glob import glob

lenlimit=10000
outdir=argv[1]
allprots=glob("%s/fasta/*"%outdir)

def failWrite(prot):
    name=prot.split("/")[-1]
    fasta=open(prot+"/"+name+".fasta")
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

for prot in allprots:
    target=prot.split("/")[-2]    
    try:
        f=open(prot+"fu.txt")
        fures=f.readline()
        if fures[:2]=="No":
            print("Target %s's FUpred failed"%target)
            failWrite(prot)
        f.close()
    except:
        print("Target %s has no FUpred result"%target)
        failWrite(prot)
        continue

    domains=fures.split(";")[:-1]
    name=prot.split("/")[-1]
    f=open(prot+"/"+name+".fasta")
    f.readline()
    seq=""
    for line in f:
        seq+=line.strip()
    for i in range(len(domains)):
        g=open(prot+"%d.fasta"%(i+1),"w")
        g.write(">%s|%d\n"%(target,i+1))
        domainseq=""
        for part in domains[i].split(","):
            bounds=[int(j) for j in part.split("-")]
            domainseq+=seq[bounds[0]-1:bounds[1]]
        if len(domainseq) > lenlimit:
            print("Warning: Domain %s_%d is longer than PEPPI can handle.  Only the first %d amino acids were written."%(target,i+1,lenlimit))
        g.write(fill(domainseq[:lenlimit]))
        g.close()



