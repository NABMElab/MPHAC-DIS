
file = [2 5 10];
ratio_change = [1:1:100]

for i = 1:length(file)
    for j = 1:length(ratio_change)
        coverage1 = 100;
        coverage2 = 100*ratio_change(j)*file(i)
        coverage3(j,i) = (coverage2+coverage1)/file(i)
    end
end