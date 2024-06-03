clear all
addpath 'C:\Users\94249\Desktop\NABME'
Yield = 0.01:0.01:0.99;
cycle=20;
conc = 1e-14;
conc2 = 4e-7;
Temp = 63;


for i = 1:length(Yield)
    i
    dG(i) = Yield2dG(Yield(i),conc, conc2, Temp);
    eff1(i) = onecycle_extension(dG(i))-1;
end

for i = 1:length(Yield)
    for j = 1:length(Yield)
        eff3(i,j) = equaleff(eff1(i),eff1(j),cycle);
    end
end

heatmap(eff3)
grid off
colormap    jet

writematrix(eff3,'Yiled-eff.csv')




