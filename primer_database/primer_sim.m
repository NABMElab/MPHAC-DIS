parpool('Processes',11);

clear all

seqlen = [15:1:30];
seqnum = 100000;

for i = 1:seqnum
    if mod(i,100) == 0
        i
    end
    templen = randsample(seqlen,1);
    data_ori{i,1} = cell2mat(randseq_fixlength(templen,1,0,1));
end


%dG check

data_ori = PrimerFeature_dg(data_ori);

dg = cell2mat(data_ori(:, 3));
rowsToKeep_dg = dg >= -12.5 & dg <= -10.5;

data_dg = data_ori(rowsToKeep_dg,:);

%GC check
for i = 1:length(data_dg)
    
    data_dg{i,4} = GCcount(data_dg{i,1});
    data_dg{i,5} = homopolymer(data_dg{i,1});
end

gc = cell2mat(data_dg(:, 4));
rowsToKeep_gc = gc >= 0.45 & gc <= 0.55;
data_gc = data_dg(rowsToKeep_gc,:);

%homo check
homo = cell2mat(data_gc(:, 5));

rowsToKeep_homo = homo == 0;

data_homo = data_gc(rowsToKeep_homo,:);

%dimer check
primer = data_homo(:,1);

badness = MPrimerCheck_parallel(primer,'-LA',3);

for i = 1:length(data_homo)
    i
    data_homo{i,6} = mean(badness(:,i));
    data_homo{i,7} = max(badness(:,i));
end

dimer = cell2mat(data_homo(:, 7));

rowsToKeep_dimer = dimer <= 1;

data_dimer= data_homo(rowsToKeep_dimer,:);


delete(gcp('nocreate'));



