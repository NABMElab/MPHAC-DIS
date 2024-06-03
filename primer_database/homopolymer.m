function quality = homopolymer(seq)
    if (~contains(seq,'AAAA')) && (~contains(seq,'TTTT')) && ...
        (~contains(seq,'CCC')) && (~contains(seq,'GGG'))
        quality = 0;
    else
        quality = 1;
    end
end