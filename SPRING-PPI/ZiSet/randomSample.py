from sys import argv
from random import sample

f=open(argv[1])
sampleSize=int(argv[2])
if len(argv) > 3:
    sampleCount=int(argv[3])
else:
    sampleCount=1

elements=[line.strip() for line in f]
f.close()
for i in range(sampleCount):
    g=open("sample_%d"%i,"w")
    rsamp=sample(elements,sampleSize)
    for el in rsamp:
        g.write("%s\n"%el)
    g.close()
