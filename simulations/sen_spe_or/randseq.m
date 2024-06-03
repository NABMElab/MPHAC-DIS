seqlength = 22;
pos = [1:1:seqlength]
AT = {'A','T'}
CG = {'C','G'}
pattern = {'AAAA','TTTT','CCCC','GGGG'}
count = 1;
CGcount = 0;
for CGratio = 0:0.1:1
    CGnum = round(seqlength*CGratio)
    for num = 1:100
        CGpos = randsample(pos,CGnum)
        seq{count,1} = '';
        for i = 1:seqlength
            if ismember(i,CGpos) == 1
                seq{count,1} = strcat(seq{count,1},randsample(CG,1));
            else
                seq{count,1} = strcat(seq{count,1},randsample(AT,1));
            end
        end
        count = count + 1;
    end
end
        
              
    