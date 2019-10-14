import matplotlib.pyplot as p
f=open("res.csv")
alltm=[]
rtm=[]
for line in f:
    parts=line.strip().split(",")
    at=float(parts[1])
    if at < 0:
        continue
    rt=float(parts[3])
    alltm.append(at)
    rtm.append(rt)

f.close()
p.hist(alltm)
p.show()
p.hist(rtm)
p.show()
