addpath('C:\Users\94249\Desktop\NABME')
seq = readcell('80plex.csv')
count = 1
for i = 1:length(seq)
    seq{count,2} = ToeholdDG(seq{count,1},'-T',63,'-S',0.22,'-DD');
    seq{count,3} = onecycle_extension(seq{count,2});
    %seq{count,4} = fullaccesss(10000,2,seq{count,3},20,100)
    %seq{count,4} = ceil(count/10)
    count = count + 1
end

for i = 1:80
    eff1 = seq{i,3};
    eff2 = seq{i+80,3};
    eff = equaleff(eff1,eff2,20);
    seq{i,5} = eff;
    seq{i+80,5} = eff;
    eff1 = seq{i+160,3};
    eff2 = seq{i+160+80,3};
    eff = equaleff(eff1,eff2,20);
    seq{i+160,5} = eff;
    seq{i+160+80,5} = eff;
end