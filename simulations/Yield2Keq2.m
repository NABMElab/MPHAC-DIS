function [Keq] = Yield2Keq2(yield,conc, varargin)
% [dG] = Keq2dG(Keq, <Temp>)
%Calculates dG from input Keq, and optional temperature
%Temperature defaults to 37 C.
Temp = 37+273.15;
if (length(varargin) == 1)
Temp = varargin{1} + 273.15;
end
c2 = conc*yield;
c1 = conc*(1-yield);
Keq = c2/(c1*c1);
end