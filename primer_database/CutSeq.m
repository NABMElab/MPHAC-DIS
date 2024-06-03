function [segment] = CutSeq(seq,size)
count = 1;
for i = 1:floor(length(seq)/size)
    tempseq = seq(((i-1)*size+1):(i*size))
    if isempty(strfind(tempseq,'N'))
        segment{count,1} = tempseq;
        count = count + 1
    end
end