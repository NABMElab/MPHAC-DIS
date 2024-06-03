function [yield] = CalcYield2(dG, Temp, conc, conc2)
Keq = dG2Keq(dG, Temp);
if min(roots([Keq-1, -Keq*(conc+conc2), Keq*conc*conc2]))>0
yield = min(roots([Keq-1, -Keq*(conc+conc2), Keq*conc*conc2]))/min(conc, conc2);
else
yield = max(roots([Keq-1, -Keq*(conc+conc2), Keq*conc*conc2]))/min(conc, conc2);
end
end
% calculate yield of strand displacement (A + BC → AB + C)