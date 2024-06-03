addpath('/data/wengz/DNA_storage/matlab_script');
fid = fopen('/data/wengz/DNA_storage/primer_oligo_sequence/ref/coverage_ref.txt', 'r');
tempstr = fgetl(fid);
a = 1
while ischar(tempstr)
    ref{a} = tempstr(3:end);
    tempstr = fgetl(fid);
    a = a + 1;
end
fclose(fid)

 DDJ_FP = 'AAAGACGTCACAGCAAGGTTCAAATCATTCTCTCCTATCTCATC';% FP+b from second base
 DDJ_RP = 'TGAAGTGATGTGACAGCTCAGC';
 SZJ_FP = 'TAGGAGAGATTGGGCTAGAGAGATAATTGAGTGTCATCAGAACTAGAT';% FP+b from second base
 SZJ_RP = 'ATACCATGAGAGGCAGAGTGTGG';
 SZBF_FP = 'CTTTATCAGACACAGTTATGTGCTGGAAAGAGCATAAATTTTGGAAT';% FP+b from second base
 SZBF_RP = 'ACTGTGTTCTGTCACCTCTGTG';
 LY_FP = 'ATGGGACTCCAATGCAAAACTCAATGTATCAGTGTGAGGATGT';% FP+b from second base
 LY_RP = 'TGATGTGGCATATCCTTCAGTAACTT';
 AL_FP = 'ATGGGACTCCAATGCAAAACTCAATGTATCAGTGTGAGGATGT';% FP+b from second base
 AL_RP = 'AGGACACTTTGAGATCTGGCTTC';
 ML_FP = 'CAATGGGAGTCACTGCTGCAGGCGCCCTGT';% FP+b from second base
 ML_RP = 'CGACGCTGATCCAGTGTACG';
 SJTU_FP = 'TCTTCCTCTCACATCTTTATTTAACCCATTAGAAAATCCTATCAGCTCT';% FP+b from second base
 SJTU_RP = 'AGGCACAAGCTGGCAGTG';
 IHAD_FP = 'TCATCTGTAAAGCAGGGAGAGAACCTCCTCCCTCACAGA';% FP+b from second base
 IHAD_RP = 'CTGCTGGATATCTGATGGCTGT';


fid = fopen(libname, 'r');
tempstr = fgetl(fid);
count = zeros(1,length(ref));
seq = 0;

while ischar(tempstr)
    tempseq = fgetl(fid);
    fgetl(fid); fgetl(fid);
    tempstr = fgetl(fid);
    seq = seq + 1
    if (contains(tempseq(1:70),ML_FP)) & (contains(tempseq(90:end),ML_RP))
        pos = strfind(tempseq(1:70),ML_FP);
        [~,index] = ismember(tempseq(pos:pos+119),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),LY_FP)) & (contains(tempseq(90:end),LY_RP))
        pos = strfind(tempseq(1:70),LY_FP);
        [~,index] = ismember(tempseq(pos:pos+138),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),AL_FP)) & (contains(tempseq(90:end),AL_RP))
        pos = strfind(tempseq(1:70),AL_FP);
        [~,index] = ismember(tempseq(pos:pos+135),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end 
     elseif (contains(tempseq(1:70),SJTU_FP)) & (contains(tempseq(90:end),SJTU_RP))
        pos = strfind(tempseq(1:70),SJTU_FP);
        [~,index] = ismember(tempseq(pos:pos+136),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),SZBF_FP)) & (contains(tempseq(90:end),SZBF_RP))
        pos = strfind(tempseq(1:70),SZBF_FP);
        [~,index] = ismember(tempseq(pos:pos+138),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end  
  
    elseif (contains(tempseq(1:70),DDJ_FP)) & (contains(tempseq(90:end),DDJ_RP))
        pos = strfind(tempseq(1:70),DDJ_FP);
        [~,index] = ismember(tempseq(pos:pos+135),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),IHAD_FP)) & (contains(tempseq(90:end),IHAD_RP))
        pos = strfind(tempseq(1:70),IHAD_FP);
        [~,index] = ismember(tempseq(pos:pos+130),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),SZJ_FP)) & (contains(tempseq(90:end),SZJ_RP))
        pos = strfind(tempseq(1:70),SZJ_FP);
        [~,index] = ismember(tempseq(pos:pos+140),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end  
 
    end

end
fclose(fid)
count
writematrix(count,strcat('number',num2str(libnumber),'.csv'))
coverage = tabulate(count)
writematrix(coverage,strcat('coverage',num2str(libnumber),'.csv'))
