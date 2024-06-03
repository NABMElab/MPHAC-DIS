function [Tm] = FindTm(seq, concentration, varargin)
if (length(varargin) == 1)
concentration2 = varargin{1};
else
concentration2 = concentration;
end
R = 8.314/4180;
Tm = DuplexDH(seq) / (DuplexDS(seq)/1000 + ...
R*log(concentration2 - concentration/2)) - 273.15;
end
