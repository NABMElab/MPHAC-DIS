function primer = FindFP(sequence,energy,temperature)
%function for finding forward primer at specific location of a given
%sequence. 

%addpath('/Users/jiaming/Dropbox (Nablab)/MatlabCodes/PrimerDesign/EnergyCalculations');
%addpath('/Users/jiaming/Dropbox (Nablab)/MatlabCodes/PrimerDesign/MATLAB_Codes');

% ind2 = length(sequence)-20;&& ind1 = ind2-16; for outer primer design

% ind2 = length(sequence)-6;&& ind1 = ind2-16; for inner primer design

ind2 = length(sequence)-20;
ind1 = ind2-16;
seq = sequence(ind1:ind2);
dg = ToeholdDG(seq,'-T',60,'-S',0.18,'-DD');
count=0;
while dg>-energy(1)||dg<-energy(2)
    count = count+1;
    if count>=7
        ind2 = ind2+1;
        count = 0;
    end
    if dg >-energy(1)
        ind1 = ind1-1;
    else
        ind1 = ind1+2;
    end
    seq = sequence(ind1:ind2);
    dg = ToeholdDG(seq,'-T',60,'-S',0.18,'-DD');
end
primer{1,1} = seq;
primer{2,1} = [ind1 ind2];
primer{3,1} = dg;
end


