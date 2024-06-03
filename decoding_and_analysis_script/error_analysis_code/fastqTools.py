# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:20:03 2023

@author: AA
"""
import pandas as pd
import csv

def readFastq(fastqFile, compression='gzip'):
    o = pd.read_csv(fastqFile, header=None, sep='\t', compression=compression)
    df = pd.DataFrame(o.values.reshape(-1,4))
    df.columns = ['info','seq','addi','quality']
    print('[^~^] FASTQ transformed.')
    return(df)

def writeFastq(out, fastqOut):
    w = pd.DataFrame(out.values.reshape(-1,1))
    w.to_csv(fastqOut, header=0, index=0, sep='\t', quoting=csv.QUOTE_NONE)
    print('[^~^] FASTQ wrote.')

def countSeq(df, seq):
    count = df[df['seq'].str.contains(seq)]
    print("The seq: "+str(count.shape[0]))
    return(count.shape[0])

def truncate_text(row):
    length = row['seqlength']
    truncated_text = row['quality'][:length]
    return truncated_text

def trimSeq(df, adapter, part=0):
    trimmed = df.copy(deep=False)
    trimmed['seq'] = trimmed['seq'].str.split(adapter, expand=True)[part]
    trimmed['seqlength'] =  trimmed['seq'].str.len()
    trimmed['quality'] = trimmed.apply(truncate_text, axis=1)
    trimmed = trimmed.iloc[:,0:4]
    print('[^~^] FASTQ trimmed.')
    return(trimmed)

def reverseComplement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'T', 't':'A', 'c':'G', 'g':'C','N':'N','n':'n','-':'-','.':'.', '?':'?','_':'_'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases[::-1])

def mutClassify(base):
    mut = {'AT':'A>T', 'AC':'A>C', 'AG':'A>G', 'AN':'A>N', 'TA':'T>A', 'TC':'T>C', 'TG':'T>G', 'TN':'T>N',
           'CT':'C>T', 'CA':'C>A', 'CG':'C>G', 'CN':'C>N', 'GA':'G>A', 'GC':'G>C', 'GT':'G>T', 'GN':'G>N'}
    return(mut[base])