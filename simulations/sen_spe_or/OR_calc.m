clear all
dG = 0:-0.1:-16;
cycle=20;
ddG = 4;

dG2 = dG + ddG;

for i = 1:length(dG)
        eff1(i) = onecycle_extension(dG(i))-1;
end

for i = 1:length(dG2)
        eff2(i) = onecycle_extension(dG2(i))-1;
end



for i = 1:length(dG)
    for j = 1:length(dG)
        eff3(i,j) = equaleff(eff1(i),eff1(j),cycle);
    end
end

for i = 1:length(dG2)
    for j = 1:length(dG2)
        eff4(i,j) = equaleff(eff2(i),eff2(j),cycle);
    end
end

%spe = eff3./eff4;

%spe = spe./(max(max(spe))+1)
spe = eff3./(eff3+eff4);
sen = eff3;

sen2 = 1-sen;
spe2 = 1-spe;
OR = (sen./sen2)./(spe2./spe)


writematrix(OR,'or-4.csv')
writematrix(spe,'spe-4.csv')
%writematrix(eff3,'eff.csv')

%%

heatmap(OR)
grid off
colormap    jet
h=gca
YourYticklabel=cell(size(h.YDisplayLabels))
[YourYticklabel{:}]=deal('');
h.YDisplayLabels=YourYticklabel
h.XDisplayLabels=YourYticklabel
saveas(gcf, 'OR-DDG8.jpg')

heatmap(spe)
grid off
colormap    jet
h=gca
YourYticklabel=cell(size(h.YDisplayLabels))
[YourYticklabel{:}]=deal('');
h.YDisplayLabels=YourYticklabel
h.XDisplayLabels=YourYticklabel
saveas(gcf, 'SPE-DDG8.jpg')