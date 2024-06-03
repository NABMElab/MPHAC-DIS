#!/bin/bash


foldername=20240321_coverage_FL
thresh=5
prefile=SQ24021449-library20240320-lib1-
#coverage=15

#for i in {1..8};do
for i in {3..8};do
#Êý¾ÝÎ»ÖÃ£º "/data/dataBackup/wengz/SeqData/DNA_storage/20240321_coverage_FL/SQ24021449-library20240320-lib1-1_combined_R1.fastq"
#mkdir -p /data/wengz/DNA_storage/DecodingResult/${foldername}/decoding_result_sample${coverage}/lib${i}/
#cd  /data/wengz/DNA_storage/DecodingResult/${foldername}/decoding_result_sample${coverage}/lib${i}/

mkdir -p /data/wengz/DNA_storage/DecodingResult/${foldername}/decoding_result/lib${i}/
cd /data/wengz/DNA_storage/DecodingResult/${foldername}/decoding_result/lib${i}/
cp  /data/wengz/DNA_storage/decoding_matlab/* ./

libname=/data/dataBackup/wengz/SeqData/DNA_storage/${foldername}/${prefile}${i}_combined_R1.fastq
#libname=/data/dataBackup/wengz/SeqData/DNA_storage/${foldername}/${prefile}${i}_${coverage}cov.fastq


echo ${libname}
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh}; reconstruct_SJTU; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_AL; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_ML; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_IHAD; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_LY; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_SZBF; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_SZJ; exit;"
matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};threshold=${thresh};reconstruct_DDJ; exit;"

done