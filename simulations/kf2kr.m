function [kr] = kf2kr(dG, Temp, kf)
Keq = dG2Keq(dG, Temp);
kr = kf/Keq;
end