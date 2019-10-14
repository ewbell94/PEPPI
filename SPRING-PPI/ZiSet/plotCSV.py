import matplotlib.pyplot as p
from sys import argv
from scipy.stats import ttest_rel

mmpoints={}
rtmpoints={}

f=open(argv[1])
for line in f:
    parts=line.split(",")
    mmpoints[parts[0]]=[float(parts[1])]
    rtmpoints[parts[0]]=[float(parts[3])]
f.close()

f=open(argv[2])
for line in f:
    parts=line.split(",")
    mmpoints[parts[0]].append(float(parts[1]))
    rtmpoints[parts[0]].append(float(parts[3]))
f.close()

plot=False
if len(argv)>3:
    plot=True

mmdata=[]
mmdiff=[]
for key in mmpoints.keys():
    point=mmpoints[key]
    if len(point) > 1 and min(point[0],point[1]) > 0.:
        mmdata.append(point)
        mmdiff.append((key,point[1]-point[0]))
mmdiff.sort(key=lambda x: x[1])

rtmdata=[]
rtmdiff=[]
for key in rtmpoints.keys():
    point=rtmpoints[key]
    if len(point) > 1 and min(point[0],point[1]) > 0.:
        rtmdata.append(point)
        rtmdiff.append((key,point[1]-point[0]))
rtmdiff.sort(key=lambda x: x[1])

if plot:
    p.figure(figsize=[6,6])
    p.scatter([i[0] for i in mmdata],[i[1] for i in mmdata],facecolors='none',edgecolors='black')
    p.plot([i*0.1 for i in range(11)],[i*0.1 for i in range(11)],color='black',linestyle='--')
    #p.axis('equal')
    p.xlim(0.,1.)
    p.ylim(0.,1.)
    p.xlabel("TM-score No DFIRE")
    p.ylabel("TM-score wTM-score construction")
    p.savefig("mmscore.svg")
    p.show()
print(sum([i[0] for i in mmdata])/len(mmdata),sum([i[1] for i in mmdata])/len(mmdata),len(mmdata))
print(ttest_rel([i[0] for i in mmdata],[i[1] for i in mmdata]))
print(mmdiff[:10])
print(mmdiff[-10:])
print()

if plot:
    p.figure(figsize=[6,6])
    p.scatter([i[0] for i in rtmdata],[i[1] for i in rtmdata],facecolors='none',edgecolors='black')
    p.plot([i*0.1 for i in range(11)],[i*0.1 for i in range(11)],color='black',linestyle='--')
    #p.axis('equal')
    p.xlim(0.,1.)
    p.ylim(0.,1.)
    p.xlabel("rTM-score No DFIRE")
    p.ylabel("rTM-score wTM-score construction")
    p.savefig("rtmscore.svg")
    p.show()
print(sum([i[0] for i in rtmdata])/len(rtmdata),sum([i[1] for i in rtmdata])/len(rtmdata),len(rtmdata))
print(ttest_rel([i[0] for i in rtmdata],[i[1] for i in rtmdata]))
print(rtmdiff[:10])
rev=rtmdiff[-10:]
rev.reverse()
print(rev)
