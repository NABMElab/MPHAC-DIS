function [dH] = DuplexDH(seq)

Parameter = load('DnaDna.txt');
paraH = zeros(4,4);
paraS = zeros(4,4);
for i = 1:4
    for j = 1:4
        paraH(i,j) = Parameter((i-1)*4+j,2);
        paraS(i,j) = Parameter((i-1)*4+j,3);
    end
end

rep = floor(fromGCAT(seq));

dH = 0.2;
%Stack energies
for i = 1:(length(seq)-1)
    dH = dH + paraH(rep(i),rep(i+1));
end
