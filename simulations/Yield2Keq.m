function [Keq] = Yield2Keq(Yield,conc, conc2, varargin)
% [dG] = Keq2dG(Keq, <Temp>)
%Calculates dG from input Keq, and optional temperature
%Temperature defaults to 37 C.
Temp = 37+273.15;
conc = conc;
conc2 = conc2;
if (length(varargin) == 1)
Temp = varargin{1} + 273.15;
end
T1 = conc*(1-Yield);
T2 = conc2 - conc*(1-Yield);
T3 = conc*(Yield);
Keq = T3/(T1*T2);
end