


addpath('../../')
number = 100

data1 = [randseq_fixlength(15,number,0.4,0.6),randseq_fixlength(15,number,0.4,0.6)];
data2 = [randseq_fixlength(20,number,0.4,0.6),randseq_fixlength(20,number,0.4,0.6)];
data3 = [randseq_fixlength(25,number,0.4,0.6),randseq_fixlength(25,number,0.4,0.6)];
data4 = [randseq_fixenergy(number,-10.5,-12.5,0.4,0.6),randseq_fixenergy(number,-10.5,-12.5,0.4,0.6)]

data5 = [data1;data2;data3;data4];


cycle=20;
ddG = 6;

for i = 1:length(data5)
    i
    dG1(i,1) = ToeholdDG(data5{i,1},'-T',63,'-S',0.22,'-DD');
    dG2(i,1) = ToeholdDG(data5{i,2},'-T',63,'-S',0.22,'-DD');
    eff1(i,1) = onecycle_extension(dG1(i))-1;
    eff2(i,1) = onecycle_extension(dG2(i))-1;
    eff3(i,1) = equaleff(eff1(i),eff2(i),cycle);
    Eff1(i,1) = onecycle_extension(dG1(i)+ddG)-1;
    Eff2(i,1) = onecycle_extension(dG2(i)+ddG)-1;
    eff4(i,1) = equaleff(Eff1(i),Eff2(i),cycle);
end



spe = eff3./eff4;

spe = spe./(max(max(spe))+1)
spe = eff3./(eff3+eff4);
sen = eff3;

sen2 = 1-sen;
spe2 = 1-spe;
OR = (sen./sen2)./(spe2./spe)




