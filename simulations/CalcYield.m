function [yield] = CalcYield(dG, Temp, conc, conc2)
Keq = dG2Keq(dG, Temp);
yield = min(roots([Keq, -Keq*(conc+conc2)-1, Keq*conc*conc2]))/min(conc, conc2);
end

