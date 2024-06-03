#!/bin/bash

foldername=20240321_coverage_FL
prefile=SQ24021449-library20240320-lib1-


for i in {1..8};do

mkdir -p /data/wengz/DNA_storage/out/${foldername}/

libname=/data/dataBackup/wengz/SeqData/DNA_storage/${foldername}/${prefile}${i}_combined_R1.fastq
echo ${libname}

python ./error_analysis_code/step1_fuzzywuzzy_match.py /data/dataBackup/wengz/SeqData/DNA_storage/${foldername}/${prefile} ./out/${foldername}/ ${i} ${libname}
#echo 'step1 done.'
#sh ./code/step2_bwa_align.sh ./out/ ./out/
#echo 'step2 done.'
#python ./code/step3_diffComparison.py ./out/ ./out/
#echo 'step3 done.'
#"/data/dataBackup/wengz/SeqData/DNA_storage/20231102_MLtagPCR/SQ23056588-library20231101-N-lib1_combined_R1.fastq"

done