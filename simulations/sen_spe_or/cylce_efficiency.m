efficiency1 = 1;
efficiency2 = [0:0.01:1];
cycle = [0:1:40];

for i = 1:length(efficiency2)
    for j = 1:length(cycle)
        foldchange(i,j) = ((1+efficiency1)^cycle(j))/((1+efficiency2(i))^cycle(j));
    end
end

b = log10(foldchange);
heatmap(cycle,efficiency2,log10(foldchange))