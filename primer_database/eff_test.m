addpath 'C:\Users\94249\Desktop\NABME\PCR\22'
data = randseq_fixlength(20,1100,0.45,0.55)

for i = 1:length(data)
    i
    data{i,2} = ToeholdDG(data{i,1},'-T',63,'-S',0.22,'-DD');
    data{i,3} = onecycle_extension(data{i,2});
end