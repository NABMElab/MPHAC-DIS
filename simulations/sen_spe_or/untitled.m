


seq(:,1) = randseq_fixenergy(10000,-9,-10,0,1,60,0.18)
seq(:,2) = randseq_fixenergy(10000,-9,-10,0,1,70,0.18)
for i = 1:length(seq)
    seq{i,3} = ToeholdDG(seq{i,1},'-T',60,'-S',0.18,'-DD');
    seq{i,4} = ToeholdDG(seq{i,2},'-T',70,'-S',0.18,'-DD');
    seq{i,5} = length(seq{i,1});
    seq{i,6} = length(seq{i,2});
    seq{i,7} = GCcount(seq{i,1});
    seq{i,8} = GCcount(seq{i,2});
end
