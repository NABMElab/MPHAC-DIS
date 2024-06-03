function [primer_group] = PrimerFeature(primer_group)

for i = 1:length(primer_group)
    primer_group{i,3} = ToeholdDG(primer_group{i,1},'-T',63,'-S',0.22,'-DD');
end