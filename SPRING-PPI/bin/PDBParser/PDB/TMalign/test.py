import sys
sys.path.append("/nfs/amino-home/bgovi/functions/python")
sys.path.append("/nfs/amino-home/bgovi/source/PDB/")

from GetFasta import GetFasta
from pdbchains import PDBchains
from tmalign import __tmalign__

pdb1 = PDBchains()
pdb1.Read('./','PDB1.pdb', ca_only = True)
seq1 = GetFasta("PDB1.pdb")

pdb2 = PDBchains()
pdb2.Read('./', "PDB2.pdb", ca_only = True)
seq2 = GetFasta("PDB2.pdb")
x = __tmalign__(seq1,pdb1.coord[0],seq2,pdb2.coord[0])

print x
