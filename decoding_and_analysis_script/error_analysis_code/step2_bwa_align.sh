#!/bin/bash

refDir='/data/xiuxh/dnaDiagnosis/process/step2_align/ref/'  #INDEXfold
dataDir=$1
outDir=$2

for i in {1..12};do
	echo 'lib'${i}
	for lib in DDJ SZJ SZBF LY IHAD AL ML SJTU;do
	bwa mem -t 64 ${refDir}${lib} ${dataDir}lib${i}_${lib}_matched.fastq > ${outDir}MEM_lib${i}_${lib}_ori.sam
	samtools view -F 4 ${outDir}MEM_lib${i}_${lib}_ori.sam > ${outDir}MEM_lib${i}_${lib}_mapped.sam
	samtools view -f 4 ${outDir}MEM_lib${i}_${lib}_ori.sam > ${outDir}MEM_lib${i}_${lib}_unmapped.sam
	#extract FLAG + RNAME + POS + CIGAR + SEQ (col 2 3 4 6 10 of *.sam)
	cat ${outDir}MEM_lib${i}_${lib}_mapped.sam | awk -F '\t' -v OFS='\t' '{print $2,$3,$4,$6,$10}' > ${outDir}MEM_lib${i}_${lib}_mapped_clean.sam
	
	done 
done

