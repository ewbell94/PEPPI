import os
from sys import argv

os.environ["PYTHONHASHSEED"]="56709"

parts=argv[1:3]
hashdiff=str(abs(abs(hash(parts[0]))-abs(hash(parts[1]))))
hashval=int(hashdiff[-3:])
print(hashval)


