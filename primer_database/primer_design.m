parpool('Processes',6);
clear all
seq = fastaread('chr1.fasta');
segment = CutSeq(seq.Sequence,1000);
seqdata = segment(130001:230391,1);
%seqdata = segment(130001:130010,1);

amplicon = 150;
maxshift = 50;
steps = length(seqdata);

parfor i = 1:length(seqdata)
    if mod(i,100) == 0
        i
    end
    temp_P = FP_generate(seqdata{i,1},501);
    temp_P = PrimerFeature_dg(temp_P);
    [temp_P2 temp_pos temp_dg] = SelectedPrimer_dg(temp_P,1);
    temp_seq{i} = {temp_P2 temp_pos temp_dg};
end

for i = 1:length(seqdata)
    for j = 1:1
        seqdata{i,1+j} = temp_seq{i}{1}{j};
        seqdata{i,2+j} = temp_seq{i}{2}{j};
        seqdata{i,3+j} = temp_seq{i}{3}{j};
    end
end

%设计RP
temp_seq = {}
parfor i = 1:length(seqdata)
    for j = 1:1    
        if mod(i,100) == 0
        i
        end
        A = seqdata(i,:)
        temp_P = FP_generate(A{1},A{3}+amplicon);
        temp_P = PrimerFeature_dg(temp_P);
        [temp_P2 temp_pos temp_dg] = SelectedPrimer_dg(temp_P,1);
        temp_seq2{i} = {temp_P2 temp_pos temp_dg};
    end
end

for i = 1:length(seqdata)
    for j = 1:1
        seqdata{i,4+j} = temp_seq2{i}{1}{j};
        seqdata{i,5+j} = temp_seq2{i}{2}{j};
        seqdata{i,6+j} = temp_seq2{i}{3}{j};
    end
end


% 将dg转换为数值数组
dg = cell2mat(seqdata(:, 4));
% 创建一个逻辑索引，标记dg不为0的行
rowsToKeep = dg ~= 0;

data_dg = seqdata(rowsToKeep,:);
data_dg2 = data_dg(1:90000,:);


for i = 1:length(data_dg2)
    if mod(i,100) == 0
        i
    end
    data_dg2{i,5} = GCcount(data_dg2{i,2});
    data_dg2{i,6} = homopolymer(data_dg2{i,2});
end


%筛选GC含量
% 将GC含量转换为数值数组
GC = cell2mat(data_dg2(:, 5));
% 创建一个逻辑索引，标记GC含量的行
rowsToKeep_GC = GC >= 0.45 & GC <= 0.55;

data_GC = data_dg2(rowsToKeep_GC,:);


%筛选HOMOPOLYMER
% 将HOMOPOLYMER转换为数值数组
HOMO = cell2mat(data_GC(:, 6));

rowsToKeep_HOMO = HOMO == 0;

data_HOMO = data_GC(rowsToKeep_HOMO,:);

primer = data_HOMO(:,2);

badness = MPrimerCheck_parallel(primer,'-LA',3);

for i = 1:length(data_HOMO)
    i
    data_HOMO{i,7} = mean(badness(:,i));
    data_HOMO{i,8} = max(badness(:,i));
end


%筛选dimer
%将dimer转换为数值数组
dimer = cell2mat(data_HOMO(:, 8));

rowsToKeep_dimer = dimer <= 1;

data_dimer= data_HOMO(rowsToKeep_dimer,:);

delete(gcp('nocreate'));