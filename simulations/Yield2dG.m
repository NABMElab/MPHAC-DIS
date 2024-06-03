function [dG] = Yield2dG(Yield,conc, conc2, varargin)
%Calculates dG from input Keq, and optional temperature
%Temperature defaults to 37 C.
Temp = 37;
if (length(varargin) == 1)
Temp = varargin{1};
end
Keq = Yield2Keq(Yield, conc, conc2, Temp);
dG = Keq2dG(Keq, Temp);
end
