function [Keq] = dG2Keq(dG, varargin)

% [Keq] = dG2Keq(dG, <Temp>)
%Calculates Keq from input dG, and optional temperature (in Celcius)
%Temperature defaults to 37 C.

Temp = 37+273.15;

if (length(varargin) == 1)
	Temp = varargin{1} + 273.15;
end

Keq = exp(-dG * 4180 / 8.314 / Temp);

end

