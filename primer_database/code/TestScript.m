clear all
clc
% 
   primer{1} = 'AGAGGCTTGAAAAGGGCAGAG';
   primer{2} = 'TATTAGTTTGGTGACATGCATAGGACT';
 primer{3} = 'tccaatgtcagtccggcaatgtgag';
% primer{4} = 'tatcttgccgtggtcaacggactat';
% primer{5} = 'agtcactccctaagagcgatagagg';

% primer{1} = 'ccagtgttgtaggacatatattgtacc';
% primer{2} = 'ccatccccgtgtccctc';



%should be all non-capital

[~,primer,~] = xlsread('SPEtest.xlsx','NTRK1','B11:B17');

% 
% for i= 1:length(primer)
%     primer{i} = lower(primer{i});
% %     primer{i} = primer{i}(1:end-1);
% end


for i = 1:length(primer)
    primer{i} = lower(primer{i});
    primer{i} = primer{i}((primer{i} == 'g')|(primer{i} == 'c')|(primer{i} == 'a')|(primer{i} == 't'));
end

tic
badness = MPrimerCheck(primer)
toc

SummarizeBadness(primer, badness, '-top', 5);