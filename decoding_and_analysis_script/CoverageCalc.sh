#!/bin/bash


foldername=20230105_deltaG+video
prefile=R23400658-lib0105-lib
#Êý¾ÝÎ»ÖÃ£º "/data/dataBackup/wengz/SeqData/DNA_storage/20231118_diffconc/SQ23060755-library20231118-lib3-13_combined_R1.fastq"
#for i in {1..8};do
#for i in {1..2};do
for i in {6..6};do

mkdir -p /data/wengz/DNA_storage/DecodingResult/${foldername}/decoding_result/lib${i}/
cd /data/wengz/DNA_storage/DecodingResult/${foldername}/decoding_result/lib${i}/
cp  /data/wengz/DNA_storage/decoding_matlab/* ./

libname=/data/dataBackup/wengz/SeqData/DNA_storage/${foldername}/${prefile}${i}_combined_R1.fastq


echo ${libname}

matlab -nodesktop -nosplash -r "libname='${libname}'; libnumber=${i};cov_video; exit;"

done