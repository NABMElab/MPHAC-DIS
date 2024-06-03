import sys
import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import csv
from fastqTools import *

inputLoc = sys.argv[1]
outputLoc = sys.argv[2]
lib = sys.argv[3]
libname = sys.argv[4]
print(libname)
def fuzzy_match_fp(text, seq):
    text_score = {}
    for lib,p in text.items():
        score = fuzz.partial_ratio(seq, p)
        text_score[lib] = score
    text_score = sorted(text_score.items(), key=lambda x: x[1], reverse=False)
    # print(text_score)
    match = text_score[-1]
    return match[0], match[1]

def fuzzy_match_rp(text, seq):
    text_score = {}
    for lib,p in text.items():
        score = fuzz.partial_ratio(reverseComplement(seq), p)
        text_score[lib] = score
    text_score = sorted(text_score.items(), key=lambda x: x[1], reverse=False)
    # print(text_score)
    match = text_score[-1]
    trimLoc = 'No matched and deleted'
    if match[1] >= 80:
        for i in range(len(seq)):
            sub = reverseComplement(seq)[i:]
            sub_ratio = fuzz.partial_ratio(sub, rpList[match[0]])
            if sub_ratio < match[1]:
                trimLoc = i-1
                #print(trimLoc)
                break
    return match[0], match[1], trimLoc

##read fp and rp
primerFile = pd.read_csv('./code/libInfo.csv', header=0, sep=',')
primerFile = primerFile.apply(lambda x:x.astype(str).str.upper())
#primerFile['rp'] = primerFile['rp'].map(lambda x:reverseComplement(x))
#primerFile = primerFile.iloc[[0,1],:]
#print(primerFile)

fpList = primerFile.set_index(['file'])['fp'].to_dict()
rpList = primerFile.set_index(['file'])['rp'].to_dict()

##read adapter
adapterFile = pd.read_csv('./code/adapter.csv', header=0, sep=',')
adapter = reverseComplement(adapterFile.iloc[1,1])

##read fastq file
#for lib in ['lib1','lib2','lib3','lib4','lib5','lib6','lib7','lib8','lib9','lib10','lib11','lib12']:

fasDf = readFastq(libname,compression=None)
print('Original Nreads: ')
print(int(fasDf.shape[0]))
#print(fasDf['seq'].str.len().value_counts())
#trim adapter
#fasDf = trimSeq(fasDf, adapter)
#print(fasDf['seq'].str.len().value_counts())
#fasDf = fasDf.iloc[0:1000]

##20230825 modified. unique fasDf['seq'] -> new dic after fuzzy matched -> matched to ori fasDf
uniqFasDf = pd.DataFrame(fasDf['seq'].value_counts())
uniqFasDf.columns = ['count']
uniqFasDf['seq'] = uniqFasDf.index
#uniqFasDf = uniqFasDf.iloc[0:10]

uniqFasDf['matched_fp'] = uniqFasDf['seq'].map(lambda x:fuzzy_match_fp(fpList, x))
uniqFasDf['fp'] = uniqFasDf['matched_fp'].map(lambda x:x[0])
uniqFasDf['fp_score'] = uniqFasDf['matched_fp'].map(lambda x:x[1])

uniqFasDf['matched_rp'] = uniqFasDf['seq'].map(lambda x:fuzzy_match_rp(rpList, x))
uniqFasDf['rp'] = uniqFasDf['matched_rp'].map(lambda x:x[0])
uniqFasDf['rp_score'] = uniqFasDf['matched_rp'].map(lambda x:x[1])
uniqFasDf['trim_loc'] = uniqFasDf['matched_rp'].map(lambda x:x[2])

#filter
uniqFasDf['primer'] = uniqFasDf['fp'] +'_'+ uniqFasDf['rp']
uniqFasDf['matchedPrimer'] = uniqFasDf['fp'] +'_'+ uniqFasDf['rp'] +'_'+ uniqFasDf['fp_score'].map(str) +'_'+ uniqFasDf['rp_score'].map(str)+'_'+ uniqFasDf['trim_loc'].map(str)
uniqFasDf = uniqFasDf[uniqFasDf['trim_loc'] != 'No matched and deleted']
uniqFasDf_filtered = uniqFasDf.loc[uniqFasDf['primer'].isin(['AL_LY','LY_AL','AL_AL','LY_LY','DDJ_DDJ','SZBF_SZBF','SZJ_SZJ','ML_ML','SJTU_SJTU','IHAD_IHAD'])]
	
#generate dictionary
seqMatched = uniqFasDf_filtered.set_index(['seq'])['matchedPrimer'].to_dict()

##match to seq in fasDf
fasDf['matched'] = fasDf['seq'].map(seqMatched)
fasDf['info'] = fasDf['info'] +' matched:'+ fasDf['matched']

# write FASTQ file
fasDf['toFile'] = fasDf['matched'].str.split('_', expand=True)[1]
fasDf['trimLoc'] = fasDf['matched'].str.split('_', expand=True)[4]
for info in ['SZJ','DDJ','SZBF','LY','AL','ML','SJTU','IHAD']:
	needTrim = fasDf[fasDf['toFile']==info]
	if needTrim.shape[0] != 0:
		needTrim['trimLoc'] = needTrim['trimLoc'].astype(int)
		needTrim['seq'] = needTrim.apply(lambda row: row['seq'][:-row['trimLoc']] if row['trimLoc'] != 0 else row['seq'], axis=1)
		needTrim['quality'] = needTrim.apply(lambda row: row['quality'][:-row['trimLoc']] if row['trimLoc'] != 0 else row['seq'], axis=1)

		#writeFastq(fasDf[fasDf['toFile']==info][['info','seq','addi','quality']], outputLoc+lib+'_'+info+'_matched.fastq')
		writeFastq(needTrim[needTrim['toFile']==info][['info','seq','addi','quality']], outputLoc+lib+'_'+info+'_matched.fastq')