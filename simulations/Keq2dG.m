function [dG] = Keq2dG(Keq, varargin)
% [dG] = Keq2dG(Keq, <Temp>)
%Calculates dG from input Keq, and optional temperature
%Temperature defaults to 37 C.
Temp = 37+273.15;
if (length(varargin) == 1)
Temp = varargin{1} + 273.15;
end
dG = - 8.314 * Temp * log (Keq) / 4180;
end

