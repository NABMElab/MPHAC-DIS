function tempseq = randseq_fixenergy(number,maxenergy,minenergy,GC_min,GC_max,temp,salt)

ori_length = 5;
base = {'A','T','C','G'};

num = 1;
while num <= number
    tempseq{num,1} = '';
    for i = 1:ori_length
        tempseq{num,1} = [tempseq{num,1} cell2mat(randsample(base,1))];
    end

    tempenergy = ToeholdDG(tempseq{num,1},'-T',temp,'-S',salt,'-DD');

    while tempenergy > maxenergy
        tempseq{num,1} = [tempseq{num,1} cell2mat(randsample(base,1))];
        tempenergy = ToeholdDG(tempseq{num,1},'-T',temp,'-S',salt,'-DD');
    end

    if (tempenergy >= minenergy) & (GCcount(tempseq{num,1}) >= GC_min) & (GCcount(tempseq{num,1}) <= GC_max)
        num = num + 1
    end
end
