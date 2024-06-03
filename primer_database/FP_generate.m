function [FP_group] = FP_generate(seq,position)

gap = [1:15];
primer_length = [15:30];

count = 1;
for i = 1:length(gap)
    for j = 1:length(primer_length)
        FP_group{count,1} = seq((position-gap(i)-primer_length(j)+1):(position-gap(i)));
        FP_group{count,2} = (position-gap(i)-primer_length(j)+1);
        count = count + 1;
    end
end
