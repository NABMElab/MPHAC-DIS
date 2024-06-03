function [output]=GCcount(inputseq)
counta = 0;
countt = 0;
countc = 0;
countg = 0;
for i = 1:length(inputseq)
        if inputseq(i) == 'A' || inputseq(i) == 'a'
            counta = counta + 1;
        elseif inputseq(i) == 'T' || inputseq(i) ==  't'
            countt = countt + 1;
        elseif inputseq(i) == 'C' || inputseq(i) ==  'c'
            countc = countc + 1;
        elseif inputseq(i) == 'G' || inputseq(i) ==  'g'
            countg = countg + 1;
        end
end
countall = counta + countt + countc + countg;
ratioa = counta/countall;
ratiot = countt/countall;
ratioc = countc/countall;
ratiog = countg/countall;
%output = strcat('A: ',num2str(ratioa*100),'% T: ',num2str(ratiot*100),'% C: ',num2str(ratioc*100),'% G: ',num2str(ratiog*100),'%');
output = ratioc + ratiog;
    