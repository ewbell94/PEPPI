import matplotlib.pyplot as p
from numpy import std

f=open("target_list")
lens=[]
for line in f:
    try:
        g=open("%s/SPRING/TemplateSummary.txt"%line.strip())
    except:
        print("Template Summary does not exist for target %s"%line.strip())
        continue
    l=g.readlines()
    lens.append(len(l))

lens=sorted(lens)
perfects=[1 for i in lens if i==100]
print(sum(perfects),len(lens))
print(lens[:10])
print(lens[-10:])
mu=sum(lens)/len(lens)
sd=std(lens)
print(mu,sd)
#p.hist(lens)
#p.show()
