close all

addpath('/data/wengz/DNA_storage/matlab_script');
b=imread('/data/wengz/DNA_storage/primer_oligo_sequence/SJTU.png');

% 
% %b=imread('airline.bmp');
% %DDJ FP+B 46nt, 
% 
%threshold = 1;
fid = fopen(libname, 'r')
 
tempstr = fgetl(fid);
counter = 1;
counter1 = 1; 
counter2 = 1; 
counter3 = 1; 

AL_FP = 'CTCTTCCTCTCACATCTTTATTTAACC';% FP+b from second base
 AL_RP = 'AGGCACAAGCTGGCAGTG';

while ischar(tempstr)
    
    tempseq = fgetl(fid);
    
    if (strfind(tempseq(2:28),AL_FP)~=0) & (strfind(tempseq(122:end),AL_RP)~=0)
        seq{counter} = tempseq;
        counter = counter + 1;
    %if length(tempseq) == 120  && isempty(find(tempseq(46:115) == 'N', 1))
%     if length(tempseq) == 140 && isempty(find(tempseq == 'N', 1))
%         seq{counter} = tempseq;
%         counter = counter + 1;
%     elseif length(tempseq) == 140 && ~isempty(find(tempseq == 'N', 1))
%         seq1{counter1} = tempseq;
%         counter1 = counter1 + 1;
%     elseif length(tempseq) ~= 140  && isempty(find(tempseq == 'N', 1))
%         seq2{counter2} = tempseq;
%         counter2 = counter2 + 1; 
%     elseif length(tempseq) ~= 140  && ~isempty(find(tempseq == 'N', 1))
%         seq3{counter3} = tempseq;
%         counter3 = counter3+1;
     end
    
    fgetl(fid); fgetl(fid);
    tempstr = fgetl(fid);
end

fclose(fid);

%  X = [counter-1 counter1-1 counter2-1 counter3-1];
% %X = [counter-1 counter2-1+counter3-1];
% % labels = {'Exact Length no N', 'Exact Length with N', 'Different Length no N', 'Different Length with N'};
% pie(X)
% X./sum(X)*100
% totalreads = sum(X);
% 
% save('ImageSeqs.mat', 'seq');

% ------
%%
% load('ImageSeqs.mat');
% addpath('/Users/pingsong/Documents/multiple_codes/PrimerDimerCheck/Nablab_Ecalc');

num=1;

% Al, size (b,1) = hight=1=108, size(b,2) = width=2=201
for i=1:43
    for j=1:86
        Address_ref{num} = revcompseq(Remap([i,j]));
        Address_ref{num} = upper(Address_ref{num});
        num = num+1;
    end
end

% for i=1:97
%     for j=1:126
%         Address_ref{num} = revcompseq(Remap([i,j]));
%         Address_ref{num} = upper(Address_ref{num});
%         num = num+1;
%     end
% end

count=ones(1,length(Address_ref));

for i= 1:length(seq)

%         index = find(ismember(Address_ref,seq{i}(26:35)));
    [~,index] = ismember(seq{i}(52:61),Address_ref);
%     [~,index] = ismember(seq{i}(26:35),Address_ref);

    if index == 0
        continue
        
    else
        Message{index,count(index)} = seq{i}(62:121); % message to message end
        Pixel{index,count(index)} = seq{i}(52:121);
%         Message{index,count(index)} = seq{i}(36:95);

%         Message{index,count(index)} = seq{i};

        count(index) = count(index)+1;
    
    end
    
    if ~mod(i*100, length(seq))
        fprintf('.');
    end
    
end

%%

DominantMess = cell(1,14792); % 21708= 86*172
ConsensusRatio = -1*ones(1,14792);

for i=1:size(Message,1)
    AllMessPerAdd = Message(i,:);
    if ~isempty(AllMessPerAdd{1})
        MessPerAdd{i} = AllMessPerAdd(~cellfun('isempty',AllMessPerAdd));
        y = unique(MessPerAdd{i});
%   
        n = zeros(length(y), 1);
        for j = 1:length(y)
            n(j) = length(find(strcmp(y(j), AllMessPerAdd)));
        end
        
        [dominant_count, itemp] = max(n);
        [second_dominant_count, itemp2] = max(n(n~=max(n)));
%         ConsensusRatio(i) = dominant_count/length(MessPerAdd{i});
        
        if max(n) >= threshold 
            %             if length(second_dominant_count) == 0 
                DominantMess{i} = y(itemp);
                Dominant_count(i) = dominant_count;
                        
                ConsensusRatio(i) = dominant_count/length(MessPerAdd{i});

%             elseif dominant_count/second_dominant_count>=2
%                    DominantMess{i} = y(itemp);
%                    DominantMess{i} = char(DominantMess{i});
%                    SecondDominantMess{i} = y(itemp2);
%                    SecondDominantMess{i} = char(SecondDominantMess{i});
%                    Dominant_count(i) = dominant_count;
%              end
            
        end
        
    end

end

%%
%addpath('/Users/jangwonkim/Dropbox (NABLab)/MATLAB/DNA Encryption')

DominantMess22 = cell2table(DominantMess)
writetable(DominantMess22,strcat('lib',num2str(libnumber),'_SJTU_dominantmess.csv'))
a=-1*ones(86, 172 ,3);

count=1; 

for i=1:43
    for j=1:86
        
        inputseq  = char(DominantMess{count});
        count = count+1;
        
        if count > size(DominantMess,2)
            break;
        elseif ~isempty(inputseq)
            Decoded_Mess = Decode(revcompseq(inputseq),1);
            for k=1:3
                a(2*i-1,2*j-1,k) = Decoded_Mess(k);
                a(2*i,2*j-1,k) = Decoded_Mess(k+3);
                a(2*i-1,2*j,k) = Decoded_Mess(k+6);
                a(2*i,2*j,k) = Decoded_Mess(k+9);
            end
            
        end
    end
end
% 
%a2 = a(1:240, 1:320,:);
 a2 = a(1:size(b,1), 1:size(b,2),:);

a2=uint8(a2);

image(a2);
saveas(gcf,strcat('lib',num2str(libnumber),'_STJU.png'));




close all


DominantMess2 = readcell('/data/wengz/DNA_storage/primer_oligo_sequence/Original_Mess/OriginalMess_STJU.xlsx')

len = size(DominantMess2);

count = -1*ones(1,3720);

miscounter = 1;

for i = 1:len

count(i) = strcmp(DominantMess{i},DominantMess2{i});
if strcmp(DominantMess{i},DominantMess2{i}) == 0
MissmatchedMess{miscounter,1} = DominantMess2{i};
MissmatchedMess{miscounter,2} = DominantMess{i};
MissmatchedMess{miscounter,3} = num2str(i);
MissmatchedMess{miscounter,4} = Pixel{i};
miscounter = miscounter + 1;
end

end

if exist('MissmatchedMess')
MissmatchedMess2 = cell2table(MissmatchedMess)
writetable(MissmatchedMess2,strcat('SJTU_MissmatchedMess',num2str(libnumber),'.csv'))
end

dd = reshape(count,60,[]);

figure('Renderer', 'painters', 'Position', [10 10 1200 1500]);


%h = heatmap(dd,'ColorLimits',[0 1],'GridVisible','off');
h = heatmap(dd','ColorLimits',[0 1]);
h.Title = 'oligo comparison ';
saveas(gcf,strcat('lib',num2str(libnumber),'_STJU_oligocomparison.jpg'));
writematrix(dd,strcat('lib',num2str(libnumber),'_SJTU_oligocomparison.csv'))

