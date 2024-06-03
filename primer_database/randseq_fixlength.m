function seq = randseq_fixlength(length, number,GC_min,GC_max)

seqlength = length;
pos = [1:1:seqlength];
AT = {'A','T'};
CG = {'C','G'};
pattern = {'AAAA','TTTT','CCCC','GGGG'};
count = 1;
CGcount = 0;

CGnum_max = floor(seqlength*GC_max);
CGnum_min = ceil(seqlength*GC_min);
CGnum_pool = [CGnum_min:1:CGnum_max];
num = 1;
while num <= number
    CGnum = randsample(CGnum_pool,1);
    CGpos = randsample(pos,CGnum);
    seq{count,1} = '';
    for i = 1:seqlength
        if ismember(i,CGpos) == 1
            seq{count,1} = strcat(seq{count,1},randsample(CG,1));
        else
            seq{count,1} = strcat(seq{count,1},randsample(AT,1));
        end
    end
    seq{count,1};
    isempty(strfind(seq{count,1},pattern{1}));
    if isempty(strfind(cell2mat(seq{count,1}),pattern{1})) & isempty(strfind(cell2mat(seq{count,1}),pattern{2})) ...
            & isempty(strfind(cell2mat(seq{count,1}),pattern{3})) & isempty(strfind(cell2mat(seq{count,1}),pattern{4}))
        seq{count,1} = cell2mat(seq{count,1});
        num = num +1;
        count = count + 1;
    end
end
