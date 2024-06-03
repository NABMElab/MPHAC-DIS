function [ddG_mismatch] = MismatchDDG(seq1, seq2, varargin)

%Calculates Delta Delta G of a single base mismatch by calculating Delta G of seq1 pairing to complement of seq2, 
%    and then subtracting the Delta G of seq 2 pairing to complement of seq2
% varargin can include '-T (temperature)' (with temperature in Celsius).  Not specifying defaults to 37 C.
% varargin can include '-S (salinity)' (with salinity in Molar Na+).  Not specifying defaults to 1 M Na+.  
%     Don't have dependence of mismatch on salinity, so this will assume a worst case of only base stacks being affected.
% varargin can include '-RR' (for RNA/RNA), '-RD' (for RNA/DNA), and '-DD' (for DNA/DNA).  
%     Not specifying defaults to DNA/DNA.  RNA/DNA parameters not currently known.

Temp = 37+273.15;
Salinity = 1;
Type = 1;
dG_correct = 0;
dG_SNP = 0;

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

if (length(seq1) ~= 3)||(length(seq2) ~= 3)
	error('Lengths of sequences 1 and 2 must both be exactly 3nt long!\n');
end

if (upper(seq1(1)) ~= upper(seq2(1)))||(upper(seq1(2)) == upper(seq2(2)))||(upper(seq1(3)) ~= upper(seq2(3)))
	error('seq1 and seq2 must be the same in the first and last base, but differ in the 2nd!\n');
end

if (Type == 1)
	Parameter = load('DnaDnaMismatch.txt');
	dG_correct = DomainDG(seq2, '-T', Temp-273.15, '-S', Salinity);
elseif (Type == 2)
	error('RNA-DNA mismatch parameters not currently available!\n');
	%Parameter = load('RnaDnaMismatch.txt');
else
	Parameter = load('RnaRnaMismatch.txt');
	dG_correct = DomainDG(seq2, '-T', Temp-273.15, '-S', Salinity, '-RR');
end

for ii = 1:4
    for jj = 1:4
		for kk = 1:4
			paraH(ii,jj, kk) = Parameter((ii-1)*16+(jj-1)*4+kk,2);
			paraS(ii,jj, kk) = Parameter((ii-1)*16+(jj-1)*4+kk,3);
		end
    end
end

dH_SNP = paraH(floor(fromGCAT(seq1(1))), floor(fromGCAT(seq1(2))), floor(fromGCAT(revcompseq(seq2(2))))) ...
	+ paraH(floor(fromGCAT(revcompseq(seq2(3)))), floor(fromGCAT(revcompseq(seq2(2)))), floor(fromGCAT(seq1(2))));

dS_SNP = paraS(floor(fromGCAT(seq1(1))), floor(fromGCAT(seq1(2))), floor(fromGCAT(revcompseq(seq2(2))))) ...
	+ paraS(floor(fromGCAT(revcompseq(seq2(3)))), floor(fromGCAT(revcompseq(seq2(2)))), floor(fromGCAT(seq1(2))));

dG_SNP = dH_SNP - Temp * dS_SNP/1000;

ddG_mismatch = dG_SNP - dG_correct;