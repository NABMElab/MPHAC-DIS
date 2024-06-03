function [dG_toe] = ToeholdDG(Toeseq, varargin)

%Calculates Toehold energy; includes initiation terms
% toeseq is sequence of toehold plus neighbor base
% varargin can include '-T (temperature)' (with temperature in Celsius).  Not specifying defaults to 37 C.
% varargin can include '-S (salinity)' (with salinity in Molar Na+).  Not
% specifying defaults to 1 M Na+. songlab default:0.18 M.
% varargin can include '-RR' (for RNA/RNA), '-RD' (for RNA/DNA), and '-DD' (for DNA/DNA).  Not specifying defaults to DNA/DNA.
% ToeholdDG('AAAA','-T',60,'-S',0.18,'-DD')

Temp = 37+273.15;
Salinity = 1;
Type = 1;
Toelength = length(Toeseq)-1;

counter = 1;

while (length(varargin) >= counter)
	if (strcmp(varargin{counter}, '-RD'))
		Type = 2;
		counter = counter + 1;
	elseif (strcmp(varargin{counter}, '-RR'))
		Type = 3;
		counter = counter + 1;
	elseif (strcmp(varargin{counter}, '-DD'))
		Type = 1;
		counter = counter + 1;
	elseif (strcmp(varargin{counter}, '-T'))
		Temp = varargin{counter+1}+273.15;
		counter = counter + 2;
	elseif (strcmp(varargin{counter}, '-S'))
		Salinity = varargin{counter+1};
		counter = counter + 2;
	else
		error('Unexpected input token %s found!\n', varargin{counter});
	end
end

if (Type == 1)
	Parameter = load('DnaDna.txt');
elseif (Type == 2)
	Parameter = load('RnaDna.txt');
else
	Parameter = load('RnaRna.txt');
end

paraH = zeros(4,4);
paraS = zeros(4,4);

% Temp = Temp + 273.15; %Change temperature to Kelvin

for ii = 1:4
    for jj = 1:4
        paraH(ii,jj) = Parameter((ii-1)*4+jj,2);
        paraS(ii,jj) = Parameter((ii-1)*4+jj,3);
    end
end

paraS = paraS + 0.368*log(Salinity);

paraG = paraH - Temp * paraS / 1000;

Toenum = fromGCAT(lower(RNAtoDNA(Toeseq)));

%Initiation energies
if (Type == 1)
	dG_toe = 0.2 - Temp * (-5.7)/1000;
elseif (Type == 2)
	dG_toe = 1.9 - Temp * (-3.9)/1000;
else
	dG_toe = 0 - Temp * (-10.8)/1000;
end

%Stack energies
dH_toe = 0;
dS_toe = 0;

for ii = 1 : Toelength
    a = Toenum(ii);
    b = Toenum(ii+1);
    dG_toe = dG_toe + paraG(a,b);
    dH_toe = dH_toe + paraH(a,b);
    dS_toe = dS_toe + paraS(a,b);
end

%[dH_toe, dS_toe]