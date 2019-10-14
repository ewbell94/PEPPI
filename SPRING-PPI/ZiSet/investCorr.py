import matplotlib.pyplot as p
from sys import argv
from scipy.stats import pearsonr
from scipy.stats import spearmanr

ids=[]
tms=[]
scores=[]
f=open(argv[1])
datadir="/".join(argv[1].split("/")[:-1])
for line in f:
    parts=line.split(",")
    #print(parts[0])
    tms.append(float(parts[1]))
    ids.append(parts[0])
    g=open("%s/%s/SPRING/TemplateSummary.txt"%(datadir,parts[0]))
    firstline=g.readline()
    g.close()
    parts=firstline.split()
    scores.append([float(i) for i in parts[3:6]])
f.close()
weights=[1.0,12.0,1.4]
#weights=[-0.59,12.0,0.14]
scores=[[weights[0]*i[0]+weights[1]*i[1]+weights[2]*i[2]]+i for i in scores]

plot=False
if len(argv)>2:
    plot=True

for x in range(len(scores[0])):
    scorelist=[i[x] for i in scores]
    #print(pearsonr(scorelist,tms)[0])
    print(spearmanr(scorelist,tms)[0])
    if plot:
        p.scatter(scorelist,tms,facecolors='none',edgecolors='black')
        #p.plot([i*0.1 for i in range(11)],[i*0.1 for i in range(11)],color='black',linestyle='--')
        #p.axis('equal')
        #p.xlim(0.,100.)
        p.ylim(0.,1.)
        p.xlabel("score")
        p.ylabel("TM-score")
        p.show()

