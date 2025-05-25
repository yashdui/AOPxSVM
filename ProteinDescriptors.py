import pandas as pd
import numpy as np 
from collections import Counter
import os, sys, re
import itertools
import math
import random
import pickle


def read_fasta(iFasta):
    if not os.path.exists(iFasta):
        msg = 'Error: file %s does not exist.' % iFasta
        print(msg)
        return []
    with open(iFasta) as f:
        records = f.read()
        records = records.split('>')[1:]
        fasta_sequences = []
        
        for fasta in records:
            array = fasta.split('\n')
            header, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTUVWY-]', '-', ''.join(array[1:]).upper())
            header_array = header.split('|')
            name = header_array[0]
            #label = header_array[1] if len(header_array) >= 2 else '0'
            #label_train = header_array[2] if len(header_array) >= 3 else 'training'
            fasta_sequences.append([name, sequence])
        msg="read_fasta DONE!"
        print(msg)
        return fasta_sequences 
    
def toCTDC(inFasta):
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
        }
            
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
        }
    
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = ['PID']
    
    for p in property:
        for g in range(1, len(groups) + 1):
            header.append(p + '.G' + str(g)+"_CTD.C")
            
            
    fasta_list =read_fasta(inFasta)
    
    for i in fasta_list:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for p in property:
        
            c1 = Count(group1[p], sequence) / len(sequence)
            c2 = Count(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        encodings.append(code)

    encoding_array = np.array(encodings, dtype=str)
    
    CTDC_encoding=pd.DataFrame(encodings)
    
    CTDC_encoding.columns=header
    
   
    print("CTDC is DONE! CTDC_encoding=",CTDC_encoding.shape)
    
    return CTDC_encoding
    
def toCTDT(inFasta):
    group1 = {
                'hydrophobicity_PRAM900101': 'RKEDQN',
                'hydrophobicity_ARGP820101': 'QSTNGDE',
                'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
                'hydrophobicity_PONP930101': 'KPDESNQT',
                'hydrophobicity_CASG920101': 'KDEQPSRNTG',
                'hydrophobicity_ENGD860101': 'RDKENQHYP',
                'hydrophobicity_FASG890101': 'KERSQD',
                'normwaalsvolume': 'GASTPDC',
                'polarity': 'LIFWCMVY',
                'polarizability': 'GASDT',
                'charge': 'KR',
                'secondarystruct': 'EALMQKRH',
                'solventaccess': 'ALFCGIVW'
            }
    group2 = {
                'hydrophobicity_PRAM900101': 'GASTPHY',
                'hydrophobicity_ARGP820101': 'RAHCKMV',
                'hydrophobicity_ZIMJ680101': 'HMCKV',
                'hydrophobicity_PONP930101': 'GRHA',
                'hydrophobicity_CASG920101': 'AHYMLV',
                'hydrophobicity_ENGD860101': 'SGTAW',
                'hydrophobicity_FASG890101': 'NTPG',
                'normwaalsvolume': 'NVEQIL',
                'polarity': 'PATGS',
                'polarizability': 'CPNVEQIL',
                'charge': 'ANCQGHILMFPSTWYV',
                'secondarystruct': 'VIYCWFT',
                'solventaccess': 'RKQEND'
            }
    group3 = {
                'hydrophobicity_PRAM900101': 'CLVIMFW',
                'hydrophobicity_ARGP820101': 'LYPFIW',
                'hydrophobicity_ZIMJ680101': 'LPFYI',
                'hydrophobicity_PONP930101': 'YMFWLCVI',
                'hydrophobicity_CASG920101': 'FIWC',
                'hydrophobicity_ENGD860101': 'CVLIMF',
                'hydrophobicity_FASG890101': 'AYHWVMFLIC',
                'normwaalsvolume': 'MHKFRYW',
                'polarity': 'HQRKNED',
                'polarizability': 'KMHFRYW',
                'charge': 'DE',
                'secondarystruct': 'GNPSD',
                'solventaccess': 'MSPTHY'
            }

    groups = [group1, group2, group3]
    property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = ['PID']
    
    for p in property:
        for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
            header.append(p + '.' + tr)
    
    fasta_list=read_fasta(inFasta)

    for i in fasta_list:
        name, sequence= i[0], re.sub('-', '', i[1])
        code = [name]
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        
        for p in property:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                    
                if (pair[0] in group1[p] and pair[1] in group3[p]) or ( pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                    
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
                
            code = code + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
        encodings.append(code)

    #encoding_array = np.array(encodings, dtype=str)
    
    CTDT_encoding=pd.DataFrame(encodings)
    
    CTDT_encoding.columns=header
   
    print("CTDT is DONE! CTDT_encoding=",CTDT_encoding.shape)        
     
    return CTDT_encoding

def toCTDD(inFasta):

    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
        'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
        'hydrophobicity_PONP930101',
        'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
        'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = ['PID']
    for p in property:
        for g in ('1', '2', '3'):
            for d in ['0', '25', '50', '75', '100']:
                header.append(p + '.' + g + '.residue' + d)
    
    
    fasta_list=read_fasta(inFasta)

    for i in fasta_list:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for p in property:
            code = code +  Count1(group1[p], sequence) +  Count1(group2[p], sequence) +  Count1(
                group3[p], sequence)
                
        encodings.append(code)

    CTDD_encoding=pd.DataFrame(encodings)
    
    CTDD_encoding.columns=header
   
    print("CTDD is DONE! CTDD_encoding=",CTDD_encoding.shape) 
    
    return CTDD_encoding
    
def toCTD(inFasta):
    
    ctdc=toCTDC(inFasta)
    #print(ctdc.head())
    ctdt=toCTDT(inFasta)
    #print(ctdt.head())
    ctdd=toCTDD(inFasta)
    #print(ctdd.head())
    
    ctd=pd.concat([ctdc,ctdt.iloc[:,1:],ctdd.iloc[:,1:]],axis=1)
    print("CTD is DONE@ CTD_encoding=",ctd.shape)
    
    return ctd

def toBLOSUM62(inFasta):
 
    
    blosum62 = {
        'A': [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],  # A
        'R': [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
        'N': [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],  # N
        'D': [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],  # D
        'C': [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
        'Q': [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],  # Q
        'E': [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],  # E
        'G': [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],  # G
        'H': [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],  # H
        'I': [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],  # I
        'L': [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],  # L
        'K': [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],  # K
        'M': [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],  # M
        'F': [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],  # F
        'P': [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],  # P
        'S': [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],  # S
        'T': [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],  # T
        'W': [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],  # W
        'Y': [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],  # Y
        'V': [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],  # V
        '-': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # -
    }
    encodings = []
    header = ['PID']
    
    AA='ACDEFGHIKLMNPQRSTVWY-'
    
    fasta_list=read_fasta(inFasta)
    for i in range(21*20):
        header.append('blosum62.F' + str(i+1))
    
    
    fasta_list=read_fasta(inFasta)
    
    

    for i in fasta_list:
        name, sequence = i[0], i[1]
        code = [name]
        count=Counter(sequence)
        print(count)
        print(count.keys())
        for aa in AA:
            if aa in count.keys():
                code = code + [i*count[aa] for i in blosum62[aa]]
            else:
                code = code +  blosum62['-']
            
        encodings.append(code)
    
    
    blosum62_encoding=pd.DataFrame(encodings)
    print(blosum62_encoding)
    blosum62_encoding.columns=header
    print(" blosum62 is DONE!  blosum62_encoding=",blosum62_encoding.shape) 
    
    return blosum62_encoding

     
    