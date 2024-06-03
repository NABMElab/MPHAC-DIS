addpath('/data/wengz/DNA_storage/matlab_script');
fid = fopen('/data/wengz/DNA_storage/primer_oligo_sequence/ref/video_ref.txt', 'r');
tempstr = fgetl(fid);
a = 1
while ischar(tempstr)
    ref{a} = tempstr(3:end);
    tempstr = fgetl(fid);
    a = a + 1;
end
fclose(fid)


 XGHG_FP = 'TGTATATAGACGGTAAAATAAACACCAAGACGTGGTAAATATTTACCTGGT';% FP+b from third base
 XGHG_RP = 'GGAAGCGGCTAACTATGGCG';

 XGHGF_FP = 'TAACACCAGTTCTTCCTCCACTCCACCATGGCACCTATTA';% FP+b from third base
 XGHGF_RP = 'CAGTTGCTGGGCCTTAAAGG';

 DY_FP = 'AGATGCTTTAGGCTCATGAGTTAACAAGGAGATGATGTAGTGTAAAG';% FP+b from third base
 DY_RP = 'ACAAGTGACTTGCACACTCTCA';

 HORSEH_FP = 'AAGTTGATAAATTAAAGGACTAAGGCACAGAACAATCATGCAACTTGC';% FP+b from third base
 HORSEH_RP = 'TCAAATAGACCTTGCAGATCAGCT';

 HORSEM_FP = 'GTTGTGCTGTCCATTGGCTACTCAGTCTCGGCT';% FP+b from third base
 HORSEM_RP = 'TGATCCTTGTTTGGGAGACACTC';

 HORSEL_FP = 'TTTTATATGTTAGTGTCCCCATGGTATATTGTAAGTTGTAGGTACATACCC';% FP+b from third base
 HORSEL_RP = 'ACTTACAGATGCACAGCAGGAAG';


fid = fopen(libname, 'r');
tempstr = fgetl(fid);
count = zeros(1,length(ref));
seq = 0;

while ischar(tempstr)
    tempseq = fgetl(fid);
    fgetl(fid); fgetl(fid);
    tempstr = fgetl(fid);
    seq = seq + 1
    if (contains(tempseq(1:70),XGHG_FP)) & (contains(tempseq(90:end),XGHG_RP))
        pos = strfind(tempseq(1:70),XGHG_FP);
        [~,index] = ismember(tempseq(pos:pos+140),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),XGHGF_FP)) & (contains(tempseq(90:end),XGHGF_RP))
        pos = strfind(tempseq(1:70),XGHGF_FP);
        [~,index] = ismember(tempseq(pos:pos+143),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),DY_FP)) & (contains(tempseq(90:end),DY_RP))
        pos = strfind(tempseq(1:70),DY_FP);
        [~,index] = ismember(tempseq(pos:pos+138),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end 
     elseif (contains(tempseq(1:70),HORSEH_FP)) & (contains(tempseq(90:end),HORSEH_RP))
        pos = strfind(tempseq(1:70),HORSEH_FP);
        [~,index] = ismember(tempseq(pos:pos+141),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end
     elseif (contains(tempseq(1:70),HORSEM_FP)) & (contains(tempseq(90:end),HORSEM_RP))
        pos = strfind(tempseq(1:70),HORSEM_FP);
        [~,index] = ismember(tempseq(pos:pos+125),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end  
  
    elseif (contains(tempseq(1:70),HORSEL_FP)) & (contains(tempseq(90:end),HORSEL_RP))
        pos = strfind(tempseq(1:70),HORSEL_FP);
        [~,index] = ismember(tempseq(pos:pos+143),ref);
        if index == 0
        continue
        else
        count(index) = count(index)+1;
        end 
 
    end

end
fclose(fid)
count
writematrix(count,strcat('video_number',num2str(libnumber),'.csv'))
coverage = tabulate(count)
writematrix(coverage,strcat('video_coverage',num2str(libnumber),'.csv'))
