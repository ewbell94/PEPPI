#!/usr/bin/env python
"""
The interface.py module provides operation and data extraction
methods for comparing protein dimer (seperate chain) complexes.

For functions having pdbchains as input.  Refers
to a PDBchains object containing at least two chains stored.

Methods:
InterfaceContacts: Calculates residue residue contacts between two chains. 
FractionNativeContacts: Calculates fraction of native interface contacts between
    a template(model) and native structure.
GetInterface: Returns a PDBchains object containg the interface between two 
                chains.
MMalign: Complex structure aligment using TMalign.
MinTMscore: Computes the pairwise tmscore of for each chain in two complexes 
            and returns the minimum tmscore.
Many of the functions returns several items.  Many are stored and
returned as a Pandas Series which is a structured data array.  
A quick google search will inform you of how to use it.
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

import os
import numpy as np
import pandas as pd

PDB_PATH = os.path.dirname(os.path.realpath(__file__) ) + '/'
sys.path.append(PDB_PATH)
from pdbchains import PDBchains
from alignment import TMalign, NWalign, TMscore, CalculateAlignmentScores
from alignment import KabschFunction

def InterfaceContacts(pdbchainsA ,pdbchainsB, chainA = 0, chainB = 0, 
        distanceCutoff = 10, clashCutoff = 0):
    """
    This methods determines interchain CA contacts between chainA and chainB.  
    A N*M boolean matrix is returned where the N rows are the CA residues of 
    chainA and the M columns are the CA positions of chainB.  If an element is 
    true, the two residues are within the distanceCutoff.

    Arguments:
        pdbchainsA (PDBchains object): An instantiated pdbchains object 
            containing the chain(s) of interest.
        pdbchainsB (PDBchains object): An instantiated pdbchains object 
            containing the chain(s) of interest:
        chainA (int,optional): Which chain in pdbchainsA to use in the 
            calculation. Default uses the 1st chain.
        chainB (int,optional): Which chain in pdbchainsB to use in the 
            calculation. Default uses the 1st chain.
        distanceCutoff (float, optional): Distance cutoff for residues residue 
            intercahin distances to be considered in contact.
        clashCutoff (float, optional): Values less than clashCutoff are 
            considered as a clash and are removed from interface calculation.
    Returns (contactMatrix, distanceMatrix, Panda Series):
        contactMatrix (N*M numpy boolean array): The N rows of the CA residues 
            of chainA and the M columns of the CA residues of chainB.  A 
            contact is designated as True.
        distanceMatrix (N*M numpy array): All pairwise distances of the N rows 
            of the CA residues of chainA and the M columns of CA residues of 
            chainB.
        numContacts (float): Number of interface CA contacts.
        AInterfaceCa ([int]): List of chainA CA residue positions that are at 
            the interface. The list values range from 0 to the length of chainA 
        BInterfaceCa ([int]): List of chainB CA residue positions that are at 
            the interface. The list values range from 0 to the length of chainB.
    """
    #for chainA
    ca_pos   = pdbchainsA.ca_pos[chainA]
    chainAlen = len(ca_pos)
    CAcoordA = pdbchainsA.coord[chainA][ca_pos,:]
    #chainB
    ca_pos   = pdbchainsB.ca_pos[chainB]
    chainBlen = len(ca_pos)
    CAcoordB = pdbchainsB.coord[chainB][ca_pos,:]

    contact_list = []
    distance_list = []
    numberCa = CAcoordA.shape[0]
    for i in xrange(0,numberCa):
        resA = CAcoordA[i]
        #Check if resB is within distanceCutoff Angstroms from any ca residue in chainA
        squaredDistance = np.sum( (CAcoordB - resA)**2, axis = 1)
        isContact = (squaredDistance <= distanceCutoff**2) & (squaredDistance >= clashCutoff**2) 
        distance_list.append(squaredDistance**0.5)
        contact_list.append(isContact)

    distanceMatrix = np.array(distance_list)
    contactMatrix = np.array(contact_list)
    BInterfaceCa=[i for i in xrange(0,chainBlen) if True in contactMatrix[:,i] ]
    AInterfaceCa=[i for i in xrange(0,chainAlen) if True in contactMatrix[i,:] ]
    numContacts = np.sum(contactMatrix)
    contactSummary = pd.Series([numContacts,AInterfaceCa,BInterfaceCa],
                    index=["numContacts","AInterfaceCa","BInterfaceCa"]) 
    return (contactMatrix, distanceMatrix, contactSummary)

def GetInterface(pdbchainsA, pdbchainsB, chainA = 0, chainB = 0, 
        distanceCutoff = 10, clashCutoff = 0,ca_only = False):
    """ 
    This method creates a PDBchains object containing the interface between 
    chainA and chainB.

    Arguments:
        pdbchainsA (PDBchains object): An instantiated pdbchains object 
            containing the chain(s) of interest:
        pdbchainsB (PDBchains object): An instantiated pdbchains object 
            containing the chain(s) of interest:
        chainA (int,optional): The chain in pdbchainsA used to calculate the 
            interface. Default uses the first chain in the pdbchains object.
        chainB (int,optional): The chain in pdbchainsB used to calculate the 
            interface. Default uses the first chain in the pdbchains object.
        distanceCutoff (float,optional): Distance cutoff for residue residue 
            interchain distances to be considered in contact.
        clashCutoff (float, optional): Values less than clashCutoff are 
            considered as a clash and are removed from interface calculation.
        ca_only (bool): If true interface residues only consists of CA positions
    Returns:
        pdbInterface (pdbchains object): Returns a pdbchains object.  The first 
            chain contains the interface from chainA and the second chain will 
            contain the interface from chainB.  If there is no interface the 
            object will be empty. 
    """
    (contactMatrix, distanceMatrix,IntContactSummary) = InterfaceContacts(
                pdbchainsA, pdbchainsB, chainA=chainA, chainB=chainB, 
                distanceCutoff=distanceCutoff,clashCutoff=clashCutoff)
    caAInterface = IntContactSummary["AInterfaceCa"]
    caBInterface = IntContactSummary["BInterfaceCa"]
    #return empty pdbchains object if no interface between the two chains
    if not caAInterface or not caBInterface:
        return PDBchains()

    #else return a pdbchains object with the interface 
    caAresnum = []
    caBresnum = []
    for i in caAInterface:
        ca_pos = pdbchainsA.ca_pos[chainA][i]
        resnum = pdbchainsA.atom_info[chainA][ca_pos].res_num
        caAresnum.append(resnum)

    for i in caBInterface:
        ca_pos = pdbchainsB.ca_pos[chainB][i]
        resnum = pdbchainsB.atom_info[chainB][ca_pos].res_num
        caBresnum.append(resnum)

    pdbInterface = pdbchainsA.SliceResNum(chainA, caAresnum, ca_only = ca_only)
    pdbInterface.Append( pdbchainsB.SliceResNum(chainB,caBresnum,
                    ca_only = ca_only) )

    return pdbInterface

def FractionNativeContacts(pdbchainsA, pdbchainsB,chainsA=None,chainsB=None,
        chainAlignA1B1 = None, chainAlignA2B2 = None, distanceCutoff = 10, 
        clashCutoff = 0, pdbchainsAIsFile = False,  pdbchainsBIsFile = False, 
        useTMalign = False, useNWalign = False, nativeInterfaceContacts = None):
    """
    Calculates the fraction of interface distanct contacts in pdbchainsA that 
    are similar to those in pdbchainsB.  Similarity of residues is determined 
    by either sequence alignment or structure alignment of chainA1 to chainB1 
    and chainA2 to chainB2.

    Arguments:
        pdbchainsA (PDBchains object or file): Default is pdbchains object with 
            atleast two chains. if pdbchainsAIsFile = True.  pdbchainsA should 
            be the name of a file.
        pdbchainsB (PDBchains object or file): Default is pdbchains object with 
            atleast two chains. if pdbchainsBIsFile = True.  pdbchainsB should 
            be the name of a file.
        chainsA ([int,int],optional):  The two chains in pdbchainsA that will 
            be used to calculate the interface. The default argument uses the 
            first two chains stored.
        chainsB ([int,int],optional):  The two chains in pdbchainsB that will 
            be used to calculate the interface. The default argument uses the 
            first two chains stored.
        chainAlignA1B1 ([str,str],optional): The alignment of chainA1 to 
            chainB1. If no argument is given, the default setting uses NWalign 
            to align the chains.     
        chainAlignA2B2 ([str,str],optional):  The alignment of chainA2 to 
            chainB2.  If no argument is given, the default setting uses 
            Nwalign to align the chains.
        distanceCutoff (float): The interchain residue cutoff distance to be 
            considered a interface contact.
        clashCutoff (float, optional): Values less than clashCutoff are 
            considered as a clash and are removed from interface calculation.
        pdbchainsAIsFile (bool,optional): If set to true, the function expects 
            pdbchainsA to be a file and not a PDBchains object. 
        pdbchainsBIsFile (bool, optional): If set to true, the function expects 
            pdbchainsB to be a file and not a PDBchains object.
        useTMalign (bool): If set to true, TMalign is used to align the chain 
                pairs A1/B1 and A2/B2
        useNWalign (bool): If set to true, NWalign is used to align the chain 
                pairs A1/B1 and A2/B2
        nativeInterfaceContacts (): Default is None.  It can be the output from 
            InterfaceContacts.  This is used for large scale comparisions 
            against native interface so that native interface is only 
            calculated once.
    Returns:
        fractionNativeContacts (float): Ranges from 0 to 1.  This is the 
            fraction of native contacts in pdbchainsB that are also in 
            pdbchainsA.  If pdbchainsA or pdbchainsB dont have an interface 0 is
            returned. If either A1/B1 or A2/B2 dont produce an alignment 0 is 
            returned.
        accuracy (float): Ranges from 0 to 1.  This is the fraction of correct 
            predictions divided by the total number of predictions.
        f1score (float): Measure of accuracy that combines precision (accuracy) 
            and recall (fractionNativeContacts).
            f1score = 2*( (precision*accuracy)/(precision + accuracy) )
    Example: 
    """
    if chainsA is None:
        chainsA = [0,1]
    if chainsB is None:
        chainsB = [0,1]
    if pdbchainsAIsFile:
        pdbTmp = PDBchains()
        pdbchainsA = pdbTmp.Read(pdbchainsA, ca_only = True)
    if pdbchainsBIsFile:
        pdbTmp = PDBchains()
        pdbchainsB = pdbTmp.Read(pdbchainsB, ca_only = True)

    if pdbchainsA.num_chains < 2 or pdbchainsB.num_chains < 2:
        err_msg = ''.join(["Fraction of Native Contacts Requires pdbchainsA ",
                    "and pdbchainsB to have at least 2 protein chains\n"])
        sys.stderr.write(err_msg)
        sys.exit()
    #returns 0 if no interface in native or model. or no alignment found between
    #native chains and model chains.

    #check if alignment exits or if other alignment method selected.
    #if true obtain alignment between chainA1 and chainB1

    if chainAlignA1B1 is None or useTMalign or useNWalign:
        chainAlignA1B1 = ['','']
        seqA1 = pdbchainsA.sequence[chainsA[0]]
        seqB1 = pdbchainsB.sequence[chainsB[0]]
        if useTMalign and len(seqA1) > 5 and len(seqB1) >5 :
            (transformationMatrix, TMalignOut) = TMalign(pdbchainsA, pdbchainsB,
                         modelchain = chainsA[0], nativechain = chainsB[0])
            chainAlignA1B1[0] = TMalignOut['seqModel_align']
            chainAlignA1B1[1] = TMalignOut['seqNative_align']
        else:
            NWoutput = NWalign(seqA1, seqB1)
            chainAlignA1B1[0] = NWoutput['seqA_align']
            chainAlignA1B1[1] = NWoutput['seqB_align']

    #check if alignment exits or if other alignment method selected.
    #if true obtain alignment between chainA1 and chainB1
    if chainAlignA2B2 is None or useTMalign or useNWalign:
        chainAlignA2B2 = ['','']
        seqA2 = pdbchainsA.sequence[chainsA[1]]
        seqB2 = pdbchainsB.sequence[chainsB[1]]
        if useTMalign and len(seqA2) > 5 and len(seqB2) > 5:
            (transformationMatrix, TMalignOut) = TMalign(pdbchainsA, pdbchainsB,
                         modelchain = chainsA[1], nativechain = chainsB[1])
            chainAlignA2B2[0] = TMalignOut['seqModel_align']
            chainAlignA2B2[1] = TMalignOut['seqNative_align']
        else:
            NWoutput = NWalign(seqA2,seqB2)
            chainAlignA2B2[0] = NWoutput['seqA_align']
            chainAlignA2B2[1] = NWoutput['seqB_align']

    #Determine the residue positions for alignment between A1/B1 and A2/B2
    modelAlignPos  = [ [],[] ] #residue position in model where there is alignment to native.
    nativeAlignPos = [ [],[] ] #residue position in native where there is alignment to model
    chainsAlign = [ chainAlignA1B1, chainAlignA2B2]
    for chain in [0,1]:
        mcount = 0
        ncount = 0
        alignLen = len(chainsAlign[chain][0])
        for i in xrange(0, alignLen):
            if chainsAlign[chain][0][i] != '-' and chainsAlign[chain][1][i] != '-':
                modelAlignPos[chain].append(mcount)
                nativeAlignPos[chain].append(ncount)
                mcount += 1
                ncount += 1
            elif chainsAlign[chain][0][i] != '-':
                mcount += 1
            else:
                ncount += 1 
        if len(modelAlignPos[chain]) == 0:
            return 0 #no alignment between chains
    
    
    #determine the number of interface contacts in complexA (model) and complexB (native)
    (contactMatrixModel,distanceMatrix,IntModelSummary)  = InterfaceContacts(
                pdbchainsA ,pdbchainsA, chainA = chainsA[0], 
                chainB = chainsA[1], distanceCutoff = distanceCutoff,
                clashCutoff = clashCutoff)
    contactMatrixNative = None
    distanceMatrix = None
    IntNativeSummary = None
    if nativeInterfaceContacts is None:
        (contactMatrixNative,distanceMatrix,IntNativeSummary)=InterfaceContacts(
                pdbchainsB ,pdbchainsB, chainA = chainsB[0], 
                chainB = chainsB[1], distanceCutoff = distanceCutoff,
                clashCutoff = clashCutoff)
    else:
        (contactMatrixNative,distanceMatrix,IntNativeSummary) = nativeInterfaceContacts

    #count number of native contacts and predicted contacts
    nativeContactSum = IntNativeSummary["numContacts"]
    predictedContactSum = IntModelSummary["numContacts"]
    if nativeContactSum == 0 or predictedContactSum == 0:
        return 0.0, 0.0, 0.0

    #Get subset of contactMatrixModel and contactMatrixNative that have model native alignment
    #native
    nativeChain2Align = [ [i] for i in nativeAlignPos[1] ] #numpy requires list of list for submatrix slice
    nativeChain1Align = nativeAlignPos[0]
    contactMatrixNative = contactMatrixNative[ nativeChain1Align, 
                            nativeChain2Align]
    #model
    modelChain2Align = [ [i] for i in modelAlignPos[1] ] 
    modelChain1Align = modelAlignPos[0]
    contactMatrixModel = contactMatrixModel[modelChain1Align, modelChain2Align]

    #calculate model contacts that are the same as native interface contacts; then divide by total
    #of native contacts to obtain fraction of native contacts.
    nativeContacts = np.sum( np.logical_and( contactMatrixModel, contactMatrixNative ) )*1.0
    fractionNativeContacts = nativeContacts/nativeContactSum #recall
    accuracy = nativeContacts/predictedContactSum #precision, tp/tp+fp

    if accuracy == 0 or fractionNativeContacts == 0:
        return fractionNativeContacts, accuracy, 0.0
    else:
        f1score = 2*( (fractionNativeContacts*accuracy)/(fractionNativeContacts + accuracy))
        return fractionNativeContacts, accuracy, f1score

def ComplexTMscore(modelComplex, nativeComplex,applyTransform = False,
    distanceCutoff = 10):
    """
    ComplexTMscore returns several structural scores for comparing a model 
    Complex to the native Complex. All scores are calculated after a global 
    superposition of the model to the native using the TMscore algorithm
    Google yang zhang TMscore for more information.

    Arguments:
        modelComplex (PDBchains object): The PDBchains object data structure 
            from pdbchains.py.  The object requires that the model have two 
            interacting chains stored.
        nativeComplex (PDBchains object): The PDBchains object.  This is 
            the native complex.
        applyTransform (bool): If true applies TMscore transform of model 
            onto native.
        distanceCutoff (float): ca ca contact distance for interface
            consideration
    Returns transformMatrix, (panda Series):
        tmscoreNativeChain1 (float): TMscore between 1st chain of model and 
                                     native normalized by length of native chain
        tmscoreNativeChain2 (float): TMscore between 2nd chain of model and 
                                     native normalized by length of native chain
        minTMscoreNative    (float): Minimum between tmscoreNativeChain1 and 
                                     tmscoreNativeChain2
        tmscoreModelChain1  (float): TMscore between 1st chain of model and 
                                     native normalized by length of model chain
        tmscoreModelChain2  (float): TMscore between 2nd chain of model and 
                                     native normalized by length of model chain
        minTMscoreModel     (float): Minimum between tmscoreModelChain1 and 
                                     tmscoreModelChain2
        globalTMscoreNative (float): Global TMscore of model complex to native 
                                     complex normalized by native complex length
        globalTMscoreModel  (float): Global TMscore of model complex to native 
                                     complex normalized by model complex length
        rmsd1               (float): RMSD between 1st chain of model and native 
                                     in the aligned regions by TMscore Algorithm
        rmsd2               (float): RMSD between 2nd chain of model and native 
                                     in the aligned regions by TMscore Algorithm
        rmsdAligned         (float): RMSD between model and native complex in 
                                     the aligned regions by TMscore Algorithm
        rawTMscore          (float): The raw tmscore is an extension of the 
                                    TMscore.  
        rawTMscore = 2*/( (1/tmscoreNativeChain1)+(1/tmscoreNativeChain2) )

        fnat                (float): Fraction of native interface contacts
        accuracy            (float): Accuracy of interace contacts
        f1score             (float): f1score of coverage and accuracy
    """
    modelTmp = PDBchains()
    modelTmp.Append(modelComplex,chains = [0] )
    modelTmp.Append(modelComplex,chains = [1] )
    modelMon = PDBchains()
    modelMon.Append(modelComplex,chains = [0] )
    modelMon.Extend(0, modelComplex,1, useRenumber = False) 

    nativeTmp = PDBchains()
    nativeTmp.Append(nativeComplex, chains = [0])
    nativeTmp.Append(nativeComplex, chains = [1])
    nativeMon = PDBchains()
    nativeMon.Append(nativeComplex, chains = [0])
    nativeMon.Extend(0, nativeComplex,1,useRenumber = False)

    #transformMatrix, kabschOut = KabschFunction(modelMon, nativeMon)
    #rmsdGlobal = kabschOut['rmsd']

    fnat,acc,f1_score = FractionNativeContacts(modelTmp,nativeTmp)

    transformationMatrix, tmout = TMscore(modelMon, nativeMon)
    globalTMscoreModel = tmout["TMscoreModel"]
    globalTMscoreNative = tmout["TMscoreNative"]
    modelTmp.Transformation(transformationMatrix)   
    tmscoreModelChain1, tmscoreNativeChain1, rmsd1 = CalculateAlignmentScores(
                                    modelTmp, nativeTmp, 
                                    modelchain = 0, nativechain = 0)
    tmscoreModelChain2, tmscoreNativeChain2, rmsd2 = CalculateAlignmentScores(
                                    modelTmp, nativeTmp, 
                                    modelchain = 1, nativechain = 1) 
    minTMscoreNative = min(tmscoreNativeChain1, tmscoreNativeChain2)
    minTMscoreModel = min(tmscoreModelChain1, tmscoreModelChain2)

    rmsdAligned = tmout['rmsd']
    rawTMscore = 0
    if tmscoreNativeChain1 != 0 and tmscoreNativeChain2 != 0:
        rawTMscore = 2.0/((1.0/tmscoreNativeChain1)+ (1.0/tmscoreNativeChain2) )


    transformMatrix, kabschOut = KabschFunction(modelMon, nativeMon)
    rmsdGlobal = kabschOut['rmsd']

    output = pd.Series([minTMscoreNative, minTMscoreModel, globalTMscoreNative, 
                globalTMscoreModel, tmscoreNativeChain1, tmscoreNativeChain2,
                tmscoreModelChain1, tmscoreModelChain2, rmsd1, rmsd2,
                rmsdAligned,rawTMscore,fnat,acc,f1_score,rmsdGlobal], 
                index = ['minTMscoreNative', 'minTMscoreModel', 
                'globalTMscoreNative', 'globalTMscoreModel', 
                'tmscoreNativeChain1', 'tmscoreNativeChain2',
                'tmscoreModelChain1', 'tmscoreModelChain2', 
                'rmsd1', 'rmsd2','rmsdAligned','rawTMscore','fnat','acc',
                'f1score','rmsdGlobal'])

    if applyTransform:
        modelComplex.Transformation(transformationMatrix)

    return transformationMatrix, output

def InterfaceRMSD(modelComplex, nativeComplex, distanceCutoff = 10, 
        clashCutoff = 0, nativeInterface = None):
    """
    Calculate interface rmsd. The model and native complex residue numbers and 
    residue types (standard 20 aminos)  must be matching. The function removes 
    the native interface from nativeComplex.  The residues contained in 
    modelComplex that represent the native are used to form the model interface.
    The rmsd is calculated bewtween the two interfaces.

    Requirements:
        The atom numbers for both model and native must start at 1 and end at 
        the number of atoms.
    Arguments:



    Returns:
        rmsd (float): Interface Root mean square deviation

    """

    if nativeInterface is None:
        dc = distanceCutoff
        cc = clashCutoff
        nativeInterface = GetInterface(nativeComplex,nativeComplex,chainA = 0,
                            chainB = 1,distanceCutoff=dc,clashCutoff=cc)
    if nativeInterface.num_chains <2:
        return -1

    if len(nativeInterface.ca_pos[0]) < 1:
        return -1

    caAresnumNat = {}
    for i in nativeInterface.ca_pos[0]:
        res_num = nativeInterface.atom_info[0][i].res_num
        caAresnumNat[res_num] = None

    caBresnumNat = {}
    for i in nativeInterface.ca_pos[1]:
        res_num = nativeInterface.atom_info[1][i].res_num
        caBresnumNat[res_num] = None


    caAmodel = []
    caBmodel = []
    for i in modelComplex.ca_pos[0]:
        res_num = modelComplex.atom_info[0][i].res_num
        if res_num in caAresnumNat:
            caAmodel.append(res_num)

    for i in modelComplex.ca_pos[1]:
        res_num = modelComplex.atom_info[1][i].res_num
        if res_num in caBresnumNat:
            caBmodel.append(res_num)

    if len(caAmodel) < 1 or len(caBmodel) < 1:
        return -1   

    modelInterface = modelComplex.SliceResNum(0,caAmodel,ca_only = True)
    modelInterface.Append(modelComplex.SliceResNum(1,caBmodel,ca_only = True))
    

    modelMon = PDBchains()
    modelMon.Append(modelInterface,chains = [0] )
    modelMon.Extend(0, modelInterface,1, useRenumber = False)


    nativeMon = PDBchains()
    nativeMon.Append(nativeInterface, chains = [0])
    nativeMon.Extend(0, nativeInterface,1,useRenumber = False)

    transformMatrix, kabschOutput = KabschFunction(modelMon, nativeMon)
    rmsd = kabschOutput['rmsd']
    return rmsd
