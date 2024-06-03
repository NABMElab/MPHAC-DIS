

clear all

dG = [-5:0.1:5]
conc1 = 1e-7;
conc2 = 1e-7;
conc3_fold = 10.^(-1:0.05:1);
ddG = 3;
temp = 70;


for i = 1:length(dG)
    for j = 1:length(conc3_fold)

        yield1(i,j) = CalcYield3(dG(i),temp,conc1,conc2,conc2*conc3_fold(j));
        yield2(i,j) = CalcYield3(dG(i)+ddG,temp,conc1,conc2,conc2*conc3_fold(j));
        DF(i,j) = yield1(i,j)/yield2(i,j);
    end
end

heatmap(DF)
grid off
        

%%

clear all

dG = 2
conc1 = 1e-7;
conc2 = 1e-7;
conc3_fold = 2;
ddG = 3;
temp = [50:1:70];


for i = 1:length(temp)

        yield1(i) = CalcYield3(dG,temp(i),conc1,conc2,conc2*conc3_fold);
        yield2(i) = CalcYield3(dG+ddG,temp(i),conc1,conc2,conc2*conc3_fold);
        DF(i) = yield1(i)/yield2(i);
end

plot(temp,DF)



%%


clear all

dG = [-5:0.1:5]
conc1 = 1e-7;
conc2 = 1e-7;
conc3_fold = 2;
ddG = [2:1:5];
temp = 60;


for i = 1:length(dG)
    for j = 1:length(ddG)

        yield1(i,j) = CalcYield3(dG(i),temp,conc1,conc2,conc2*conc3_fold);
        yield2(i,j) = CalcYield3(dG(i)+ddG(j),temp,conc1,conc2,conc2*conc3_fold);
        DF(i,j) = yield1(i,j)/yield2(i,j);
    end
end

heatmap(DF)
grid off






