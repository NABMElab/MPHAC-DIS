function [final_seq] = Remap_script_ps(filename, varargin)

%truth_marker = 'acagggcgcctGcagcagtgactcccattggt'; % length = 32 %take
%FP+b revcompseq
%MonaLisa

truth_marker = 'AGAGCTGATAGGATTTTCTAAtGGGTTAAATAAAGATGTGAGAGGAAGAGC'; % length = 32 %take
%FP+b revcompseq
%SJTU


%truth_marker = 'gccagagaaagtaaatatgttTttgttatgacggagcagtaaattttactctg'; %
%airline

%truth_marker = 'acatcctcacactgatacattGagttttgcattggagtcccatcc'; %
%first color picture

%truth_marker = 'tttagctgtccctctaagGtacctgaatttgacctgaatttgaatgttggcagaagt';
% first moon china take



if length(varargin)>0
    rp=varargin;
else
    %rp = 'ctttcgtaaggatggttaggtttg';
    %rp = 'cgtacactggatcagcgtcg'; %?G=-12.6kcal/mol,len=20nt % MonaLisa
    %rp = 'GAAGCCAGATCTCAAAGTGTCCT'; % airline
    %rp = 'AAGTTACTGAAGGATATGCCACATCA'; %first color pic
    %rp = 'ccacactctgcctctcatggtat'; % first moon
    rp = 'CACTGCCAGCTTGTGCCT';% SJTU

end

postseq = '';

% a=imread('Bitmap images/starbucks.bmp');
 a=imread('SJTU.png');
%a=filename;

heighta = size(a, 1);
widtha = size(a, 2);

if mod(heighta, 2)
    %Height is odd, add a row of zeros;
    a(heighta+1,:,:) = zeros(widtha, 3);
    heighta = heighta+1;
end

if mod(widtha, 2)
    %Width is odd, add a column of zeros;
    a(:, widtha+1, :) = zeros(heighta, 3);
    widtha = widtha+1;
end

numseqs = 1;

for i = 1:floor(heighta/2)
    for j = 1:floor(widtha/2)
        
        temp = [a(2*i-1,2*j-1,1), a(2*i-1,2*j-1,2), a(2*i-1,2*j-1,3)];
        temp = [temp, a(2*i,2*j-1,1), a(2*i,2*j-1,2), a(2*i,2*j-1,3)];
        temp = [temp, a(2*i-1,2*j,1), a(2*i-1,2*j,2), a(2*i-1,2*j,3)];
        temp = [temp, a(2*i,2*j,1), a(2*i,2*j,2), a(2*i,2*j,3)];
        
        tempseq = Remap(temp);
        
        final_seq{numseqs} = [rp, postseq, tempseq, Remap([i, j]), truth_marker];
        
        numseqs = numseqs + 1;
        
    end
end

for i=1:(numseqs-1)
    final_seq{i} = final_seq{i};
end

%% RP in the 5' of the final order seq, so don't take revcompseq as final seq
% for i=1:(numseqs-1)
%     final_seq{i} = revcompseq(final_seq{i});
% end

end

% %  fid = fopen('bulk.csv','w'); for i=1:length(all); fprintf(fid,'%s,\n',all{i}); end;

