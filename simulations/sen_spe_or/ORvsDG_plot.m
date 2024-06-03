clear all
dG = 0:-0.1:-16;
cycle=20;
ddG = 6;

dG2 = dG + ddG;

for i = 1:length(dG)
        eff1(i) = onecycle_extension(dG(i))-1;
end

for i = 1:length(dG2)
        eff2(i) = onecycle_extension(dG2(i))-1;
end

for i = 1:length(dG)
    eff3(i) = equaleff(eff1(i),eff1(i),cycle);
end

for i = 1:length(dG2)
        eff4(i) = equaleff(eff2(i),eff2(i),cycle);
end

spe = eff3./eff4;

spe = spe./(max(max(spe))+1)
spe = eff3./(eff3+eff4);
sen = eff3;

sen2 = 1-sen;
spe2 = 1-spe;
OR = (sen./sen2)./(spe2./spe)


plot(dG,OR)
plot(dG,spe)

data222 = [dG;OR];
data333 = data222';
%writematrix(data333,'ORvsDG.csv')