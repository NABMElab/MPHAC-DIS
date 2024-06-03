function [yield] = CalcYield3(dG, Temp, conc1, conc2, conc3) % T S P
Keq = dG2Keq(dG, Temp);

if min(roots([Keq-1, -(Keq*(conc1+conc2)+(conc3-conc2)), Keq*conc1*conc2]))>0
yield = min(roots([Keq-1, -(Keq*(conc1+conc2)+(conc3-conc2)), Keq*conc1*conc2]))/min(conc1, conc2);
else
yield = max(roots([Keq-1, -(Keq*(conc1+conc2)+(conc3-conc2)), Keq*conc1*conc2]))/min(conc1, conc2);
end
end

% calculate yield of strand displacement (A + BC → AB + C)