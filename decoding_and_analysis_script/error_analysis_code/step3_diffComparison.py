import sys
import difflib
import re
import pandas as pd
import numpy as np
import os
from fastqTools import *

inputLoc = sys.argv[1]
outputLoc = sys.argv[2]

##make dictionary (lib:seq)
primerFile = pd.read_csv('libInfo.csv', header=0, sep=',')
primerFile = primerFile.apply(lambda x:x.astype(str).str.upper())
primerFile['rp'] = primerFile['rp'].map(lambda x:reverseComplement(x))
rpList = primerFile.set_index(['file'])['rp'].to_dict()

for i in ['lib1','lib2','lib3','lib4','lib5','lib6','lib7','lib8','lib9','lib10','lib11','lib12']:
#for i in ['lib4']:
    for info,rp in rpList.items():
        print('Analyzing info:')
        print(info)
    
        ##load sam file
        if os.path.getsize(inputLoc+'MEM_'+ i + '_' + info +'_mapped_clean.sam'):
        #if os.path.getsize('/data/xiuxh/dnaDiagnosis/process/step2_align/out/MEM_'+ i + '_' + info +'_mapped_clean.sam') and info !='DDJ':
            samFile = pd.read_csv(inputLoc+'MEM_'+ i + '_' + info +'_mapped_clean.sam', header=None, sep='\t')
            samDf = samFile
            samDf.columns = ['flag','rname','pos','cigar','seq']
            print('Nreads mapped to '+info)
            print(samDf.shape[0])
            print(samDf.iloc[0:5,:])
        
            ##load oligo file
            oligoFile = pd.read_csv('/data/xiuxh/dnaDiagnosis/process/step2_align/ref/oligo_ref_'+ info +'.txt', header=None, sep='\t') 
            oligoDf = pd.DataFrame(oligoFile.values.reshape(-1,2))
            oligoDf.columns = ['rname','oligoseq']
            oligoDf['rname'] = oligoDf['rname'].str.replace('>','')
            #print(oligoDf.iloc[0:5,:])
            
            ##trim adapter
            samTrimmed = samDf.copy(deep=True)
            print(samTrimmed)
            #aa = '('+rp+')'
            #samTrimmed[['c3','tmp','c4']]=samTrimmed['seq'].str.split(aa,1,expand=True)
            #samTrimmed['tmp']=samTrimmed['tmp'].fillna('')
            #samTrimmed['seq']=samTrimmed.pop('c3') + samTrimmed.pop('tmp')
            #samTrimmed=samTrimmed.drop('c4',axis=1)
        
            ##count length of reads & output the distribution 
            bpsDf = samTrimmed.copy(deep=True)
            bpsDf['seqLength'] = bpsDf['seq'].str.len() 
            bpsCount = bpsDf['seqLength'].value_counts()
            mostBp = bpsCount.index[0]
            OutBps = pd.DataFrame(bpsCount)
            OutBps.to_csv(outputLoc+'seqLength_'+i+'_'+info+'.csv')
        
            matchedDf = pd.merge(samTrimmed, oligoDf, how='left', on='rname')
            print(matchedDf.columns.values)
            matchedDf1 = matchedDf['seq'].tolist()
            matchedDf2 = matchedDf['oligoseq'].tolist()
            matchedDf3 = []
        
            for r in range(matchedDf.shape[0]):
                #print(r)
                a = matchedDf1[r]
                b = matchedDf2[r]
                #s = difflib.SequenceMatcher(None, a, b, autojunk=False)
                s = difflib.SequenceMatcher(None, b, a, autojunk=False)
                diffOri = s.get_opcodes()
                diffTra = list(np.ravel(diffOri))
                diffFin = ' '.join(diffTra)  
                matchedDf3.append(diffFin)            

            print(len(matchedDf2))
            print(len(matchedDf3))   
            #matchedDf3 = pd.DataFrame({'seq':matchedDf1,'oligoseq':matchedDf2,'result':matchedDf3})
                    
            outFile = open(outputLoc+'diffComparison_'+i+'_'+info+'.csv', 'w')
            outFile.write('ngsSeq'+','+'oligoSeq'+','+'comparison'+'\n')
            for a in range(0,len(matchedDf3),1):
                outFile.write(matchedDf1[a]+','+matchedDf2[a]+','+matchedDf3[a]+'\n')
            outFile.close()
