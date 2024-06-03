parpool('Processes',12);
clear all
seq = fastaread('chr1.fasta');
segment = CutSeq(seq.Sequence,1000);
seqdata = segment(130001:230391,1);
seqdata = segment(130001:130002,1);

amplicon = 150;
maxshift = 50;
steps = length(seqdata);

tic
parfor i = 1:length(seqdata)
    i
    temp_P = FP_generate(seqdata{i,1},501);
    temp_P = PrimerFeature(temp_P);
    [temp_P2 temp_pos] = SelectedPrimer(temp_P,1);
    temp_seq{i} = {temp_P2 temp_pos};
end

num = 1;

for i = 1:length(seqdata)
    for j = 1:1
        seqdata{i,1+j} = temp_seq{i}{1}{j};
        seqdata{i,2+j} = temp_seq{i}{2}{j};
        [tempa tempb] = generate_20nt(seqdata{i,1},seqdata{i,2+j},maxshift);
        seqdata{i,5+j} = tempa;
        seqdata{i,6+j} = tempb;
    end
end

temp_seq = {}
parfor i = 1:length(seqdata)
    for j = 1:1    
        i
        A = seqdata(i,:)
        temp_P = FP_generate(A{1},A{3}+amplicon);
        temp_P = PrimerFeature(temp_P);
        [temp_P2 temp_pos] = SelectedPrimer(temp_P,1);
        temp_seq{i} = {temp_P2 temp_pos};
    end
end



for i = 1:length(seqdata)
    for j = 1:1
        seqdata{i,3+j} = temp_seq{i}{1}{j};
        seqdata{i,4+j} = temp_seq{i}{2}{j};
        [tempa tempb] = generate_20nt(seqdata{i,1},seqdata{i,4+j},maxshift);
        seqdata{i,7+j} = tempa;
        seqdata{i,8+j} = tempb;
    end
end

toc
delete(gcp('nocreate'));