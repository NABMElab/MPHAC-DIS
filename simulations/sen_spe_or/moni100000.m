addpath('C:\Users\94249\Desktop\NABME')
seq = data
count = 1
for i = 1:length(seq)
    seq{count,3} = ToeholdDG(cell2mat(seq{count,1}),'-T',63,'-S',0.22,'-DD');
    seq{count,4} = ToeholdDG(cell2mat(seq{count,2}),'-T',63,'-S',0.22,'-DD');
    seq{count,5} = onecycle_extension(seq{count,3});
    seq{count,6} = onecycle_extension(seq{count,4});
    seq{count,7} = equaleff(seq{count,5},seq{count,6},20);
    %seq{count,4} = fullaccesss(10000,2,seq{count,3},20,100)
    %seq{count,4} = ceil(count/10)
    count = count + 1
end

count = 1
for i = 1:length(seq)
    seq{count,1} = cell2mat(seq{count,1});
    seq{count,2} = cell2mat(seq{count,2});
    count = count + 1
end

