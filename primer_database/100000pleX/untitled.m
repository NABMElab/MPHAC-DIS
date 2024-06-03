parpool('Processes',36);

clear all

addpath('/data/wengz/NABME/');
addpath('/data/wengz/NABME/primer/primer_database/code/');
addpath(' /data/wengz/NABME/primer/primer_database/')
addpath(' /data/wengz/NABME/PCR/22/')
seqlen = 20;
seqnum = 4000000;
parfor i = 1:seqnum
    if mod(i,10000) == 0
        i
    end
    tempseq = randseq_fixlength2(seqlen,1,0,1);
    data_ori{i,1} = tempseq{1}{1};
end


%GC check
for i = 1:length(data_ori)
    
    data_ori{i,4} = GCcount(data_ori{i,1});
    data_ori{i,5} = homopolymer(data_ori{i,1});
end

gc = cell2mat(data_ori(:, 4));
rowsToKeep_gc = gc >= 0.45 & gc <= 0.55;
data_gc = data_ori(rowsToKeep_gc,:);

%homo check
homo = cell2mat(data_gc(:, 5));

rowsToKeep_homo = homo == 0;

data_homo = data_gc(rowsToKeep_homo,:);

data_final1 = data_homo(1:100000,:);
data_final2 = data_homo(100001:200000,:);
data_final = [data_final1 data_final2]



parfor i = 1:length(data_final)
    if mod(i,1000) == 0
        i
    end
    data_final{i,11}= ToeholdDG(data_final{i,1},'-T',63,'-S',0.22,'-DD');
    data_final{i,12}= ToeholdDG(data_final{i,6},'-T',63,'-S',0.22,'-DD');
    data_final{i,13}= onecycle_extension(data_final{i,11});
    data_final{i,14}= onecycle_extension(data_final{i,12});
    data_final{i,15}= equaleff(data_final{i,13},data_final{i,14});
end


parfor i = 1:length(data_final)
    if mod(i,1000) == 0
        i
    end
    temp_result(i,:) = calc_eff(data_final(i,:));
end

data_final3 = [data_final temp_result];
writecell(data_final3,'data_final_FL.csv');
delete(gcp('nocreate'));