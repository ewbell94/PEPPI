import matplotlib.pyplot as p
import numpy as np
from scipy.optimize import minimize
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from os.path import isfile
from math import sqrt

scores=[]
tms=[]
spearmans=[]

def calcSpearman(weights,scores,tms):
    springscores=[]
    print(weights)
    for scoreset in scores:
        springscores.append(weights[0]*scoreset[0]+weights[1]*scoreset[1]+weights[2]*scoreset[2])
    spearmans.append(spearmanr(springscores,tms)[0])
    return -1*spearmans[-1]

def calcPearson(weights,scores,tms):
    springscores=[]
    for scoreset in scores:
        springscores.append(weights[0]*scoreset[0]+weights[1]*scoreset[1]+weights[2]*scoreset[2])
    return pearsonr(springscores,tms)[0]

def plotCorrelation(weights,scores,tms):
    springscores=[]
    for scoreset in scores:
        springscores.append(weights[0]*scoreset[0]+weights[1]*scoreset[1]+weights[2]*scoreset[2])
    p.scatter(springscores,tms,facecolor='none',edgecolor='black')
    p.show()

f=open("target_list")
for target in f:
    target=target.strip()
    if not isfile("%s/SPRING/TemplateSummary.txt"%target) or not isfile("%s/SPRING/tms.txt"%target):
        continue
    print(target)
    a=open("%s/SPRING/TemplateSummary.txt"%target)
    b=open("%s/SPRING/tms.txt"%target)
    for aline in a:
        bline=b.readline()
        aparts=aline.split()
        bparts=bline.split(",")
        if bparts[0].split('-')[0]!=aparts[0] or bparts[0].split('-')[1]!=aparts[1]:
            print("AAAAAAA")
            continue
        tms.append(float(bparts[1]))
        scores.append([float(aparts[3]),float(aparts[4]),float(aparts[5])])
    a.close()
    b.close()

'''
p.scatter([i[0] for i in scores],tms)
print(spearmanr([i[0] for i in scores],tms)[0])
p.show()
p.scatter([i[1] for i in scores],tms)
print(spearmanr([i[1] for i in scores],tms)[0])
p.show()
p.scatter([i[2] for i in scores],tms)
print(spearmanr([i[2] for i in scores],tms)[0])
p.show()
'''

w=[1.,12.,1.4]
result=minimize(calcSpearman,w,args=(scores,tms),method='Powell')
#result=minimize(calcSpearman,w,args=(scores,tms))
if result.success:
    print("Successfully optimized")
    print(result.x)
    print(-1*calcSpearman(result.x,scores,tms))
    print(calcPearson(result.x,scores,tms))
    #plotCorrelation(result.x,scores,tms)
else:
    print("Weights not minimized")
#p.plot([i for i in range(len(spearmans))],spearmans)
#p.show()

'''
smat=np.asmatrix([i+[1.] for i in scores])
tmat=np.asmatrix(tms)
w=np.linalg.pinv(smat)*tmat.T
print(w)
print(-1*calcSpearman([w[0,0],w[1,0],w[2,0]],scores,tms))
print(calcPearson([w[0,0],w[1,0],w[2,0]],scores,tms))
'''
