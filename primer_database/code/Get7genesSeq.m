close all
clear all

addpath '/Users/pingsong/Dropbox (NABLab)/NABLab/MATLAB/Nablab_Ecalc'

[~,FGFR1] = xlsread('20201020_7genesSeq.xlsx','FGFR1','B1:B54'); %-
[~,FGFR2] = xlsread('20201020_7genesSeq.xlsx','FGFR2','B1:B54'); %-
[~,FGFR3] = xlsread('20201020_7genesSeq.xlsx','FGFR3','B1:B54'); %+
[~,FGFR4] = xlsread('20201020_7genesSeq.xlsx','FGFR4','B1:B54'); %+

[~,NTRK1] = xlsread('20201020_7genesSeq.xlsx','NTRK1','B1:B48'); %+
[~,NTRK2] = xlsread('20201020_7genesSeq.xlsx','NTRK2','B1:B57'); %+
[~,NTRK3] = xlsread('20201020_7genesSeq.xlsx','NTRK3','B1:B57'); %-

Seq_1 = FGFR1(2:3:54);
Seq_2 = FGFR2(2:3:54);
Seq_FGFR3 = FGFR3(2:3:54);
Seq_FGFR4 = FGFR4(2:3:54);

Seq_NTRK1 = NTRK1(2:3:48);
Seq_NTRK2 = NTRK2(2:3:57);
Seq_3 = NTRK3(2:3:57);

for i = 1:18
   Seq_FGFR1{i} = revcompseq(Seq_1{i});
   Seq_FGFR2{i} = revcompseq(Seq_2{i});
   
end

Seq_FGFR1 = Seq_FGFR1';
Seq_FGFR2 = Seq_FGFR2';


for j = 1:19
    Seq_NTRK3{j} = revcompseq(Seq_3{j});
end

Seq_NTRK3 = Seq_NTRK3';

FGFR = [Seq_FGFR1; Seq_FGFR2; Seq_FGFR3; Seq_FGFR4];
NTRK = [Seq_NTRK1; Seq_NTRK2; Seq_NTRK3];

for n = 1:length(FGFR)
   FGFR_revSeq{n} = revcompseq(FGFR{n});
  
end

for m = 1: length(NTRK)
    NTRK_revSeq{m} = revcompseq(NTRK{m});
end


FGFR_rev = FGFR_revSeq';
NTRK_rev = NTRK_revSeq';







