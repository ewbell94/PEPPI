# distutils: language = c++
# distutils: sources = TMalign_wrapper.cpp
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

from libcpp.string cimport string
import numpy as np

cdef extern from "TMalign_wrapper.h" namespace "TMalignC":
    cdef cppclass TMalign_wrapper:
        TMalign_wrapper() except +
        #input
        double **xa
        double **ya
        char *seqx
        char *seqy
        int xlen, ylen
        double *xresno
        double *yresno
        #output
        double rmsd_out, TM1_out, TM2_out, number_aligned
        double seqid_out
        double t0_out[3]
        double u0_out[3][3]
        string seqX_align_out, align_out, seqY_align_out
        #methods
        void allocate_memory()
        void run_tmalign()
        void char_test()    

def _tmalign(seqB, CAcoordsB, seqA, CAcoordsA):
    """
    This function is a wrapper for TMalign written in C++.
    The rotation and translation matricies retured are for
    rotation chainB onto chainA
    Arguments:
        seqA (str): Sequence of chainA
        seqB (str): Sequence of chainB 
        CAcoordsB (N*3 numpy array) Contans the CA coordinates
            of chainA
        CAcoordsB (N*3 numpy array) Contains the CA coordinates
            of chainB
    Returns:
        rotation (3*3): Rotation matrix
        translation (3*1): Translation matrix
        alignment (str): A string of ' ' and ':' representing alignment
            between chainA and chainB
        number_aligned (int): Number of residues aligned between to chains.
        rmsd (float): Root mean squared deviation of alignment.
        seqB_alignment (str): Alignment of chainB. with dashes representing
            unalinged regions.
        seqA_alignment (str): Alignment of chainA with dashes being unaligned regions.
        TM2 (float): TMscore normalized by chainB
        TM1 (float): TMscore normalized by chainA
    """
    #coordinates of the respective chains
    #B is rotated to A
    #B = X
    #A = Y
    #copies values into c module and outputs results
    cdef int i = 0
    cdef int j = 0
    cdef int lenA = len(seqA)
    cdef int lenB = len(seqB)
    transformation_matrix = np.zeros( (4,4) )
    cdef TMalign_wrapper *tm = new TMalign_wrapper()

    tm.xlen = lenB
    tm.ylen = lenA
    tm.allocate_memory()
    #copy input values for chainB
    for i in range(0,len(seqB)):
        tm.seqx[i] = ord(seqB[i])       
        tm.xresno[i] = i
        tm.xa[i][0] = CAcoordsB[i][0] 
        tm.xa[i][1] = CAcoordsB[i][1]
        tm.xa[i][2] = CAcoordsB[i][2]   
    #copy input values for chainA
    for i in range(0,len(seqA)):
        tm.seqy[i] = ord(seqA[i])
        tm.yresno[i] = i
        tm.ya[i][0] = CAcoordsA[i][0]
        tm.ya[i][1] = CAcoordsA[i][1]
        tm.ya[i][2] = CAcoordsA[i][2]
    #run tmalign.  Also deletes member pointers
    #after finsihing
    tm.run_tmalign()

    #copy output to python
    for i in range(0,3):
        transformation_matrix[i][3] = tm.t0_out[i]
        for j in range(0,3):
            transformation_matrix[i][j] = tm.u0_out[i][j]

    seqA_alignment = tm.seqY_align_out
    seqB_alignment = tm.seqX_align_out
    alignment = tm.align_out
    seqid = tm.seqid_out
    TM1 = tm.TM1_out
    TM2 = tm.TM2_out
    number_aligned = tm.number_aligned
    rmsd = tm.rmsd_out  
    #delete pointer
    del tm
    return (transformation_matrix,seqid, alignment,number_aligned,rmsd,seqB_alignment,TM2,seqA_alignment, TM1)      
