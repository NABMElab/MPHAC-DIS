close all

addpath('/data/wengz/DNA_storage/matlab_script');
% 
% %b=imread('airline.bmp');
% %DDJ FP+B 46nt, 
%
%threshold = 1; 
fid = fopen(libname, 'r')
 
tempstr = fgetl(fid);
counter = 1;

 DDJ_FP = 'ACTTTATCAGACACAGTTATGTGCT';% FP+b from second base
 DDJ_RP = 'ACTGTGTTCTGTCACCTCTGTG';
 
while ischar(tempstr)
     
     tempseq = fgetl(fid);
     
     if (strfind(tempseq(2:26),DDJ_FP)~=0) & (strfind(tempseq(120:end),DDJ_RP)~=0)
         
         seq{counter} = tempseq;
         counter = counter + 1;
 
      end
     
     fgetl(fid); fgetl(fid);
     tempstr = fgetl(fid);
 end
 
 fclose(fid);


% load('DDJ_seq.mat');
% addpath('/Users/pingsong/Documents/multiple_codes/PrimerDimerCheck/Nablab_Ecalc');

num=1;

% generate a address sequences for align, 
% size is important for the correct address position and address sequences
% i= ceil(oligos number/255), j = 255
for i=1:5
    for j=1:255
        Address_ref{num} = revcompseq(Remap([i,j]));
        Address_ref{num} = upper(Address_ref{num});
        num = num+1;
    end
end

count=ones(1,length(Address_ref));

for i= 1:length(seq)

    [~,index] = ismember(seq{i}(50:59),Address_ref);

    if index == 0
        continue
        
    else
        Message{index,count(index)} = seq{i}(60:119); % message to message end
        Pixel{index,count(index)} = seq{i}(60:119);

        count(index) = count(index)+1;
    
    end
    
    if ~mod(i*100, length(seq))
        fprintf('.');
    end
    
end

DominantMess = cell(1,1275); % 65025= 5*255
ConsensusRatio = -1*ones(1,1275);

for i=1:size(Message,1)
    
    AllMessPerAdd = Message(i,:);
    
    if ~isempty(AllMessPerAdd{1})
        
        MessPerAdd{i} = AllMessPerAdd(~cellfun('isempty',AllMessPerAdd));
        
        y = unique(MessPerAdd{i});  
        n = zeros(length(y), 1);
        
        for j = 1:length(y)
            n(j) = length(find(strcmp(y(j), AllMessPerAdd)));
        end
        
        [dominant_count, itemp] = max(n);
        [second_dominant_count, itemp2] = max(n(n~=max(n)));
        
        if max(n) >= threshold 
            
                DominantMess{i} = y(itemp);
                Dominant_count(i) = dominant_count;     
                ConsensusRatio(i) = dominant_count/length(MessPerAdd{i});
        
        end
        
    end

end

%writecell(DominantMess,strcat('lib',num2str(libnumber),'_SZBF_dominantmess.csv'))
a=-1*ones(5, 255, 2);

count=1; 

for i=1:5
    
    for j=1:255
        
        inputseq  = char(DominantMess{count});
        count = count+1;
        
        if count > size(DominantMess,2)
            break;
        elseif ~isempty(inputseq)
            
            Decoded_Mess = Decode(revcompseq(inputseq),1);
            %Decoded_Mess = Decode(inputseq,1);
            for k=1:2
                a(i,6*j-5,k) = Decoded_Mess(k);
                a(i,6*j-4,k) = Decoded_Mess(k+2);
                a(i,6*j-3,k) = Decoded_Mess(k+4);
                a(i,6*j-2,k) = Decoded_Mess(k+6);
                a(i,6*j-1,k) = Decoded_Mess(k+8);
                a(i,6*j,k) = Decoded_Mess(k+10);
            end
            
        end
    end
end

a1 = a(:,:,1);
a2 = a(:,:,2);

a3 = reshape(a1',[],1)';
a4 = reshape(a2',[],1)';

a34 = [a3' a4']';
a43 = reshape(a34,1,[]);

a43(14833:end)=[];


cc = uint8(a43);
bb = native2unicode(cc,'CP936');

writematrix(bb,strcat('lib',num2str(libnumber),'_SZBF_.dat'),'Delimiter',' ')


close all


DominantMess2 = readcell('/data/wengz/DNA_storage/primer_oligo_sequence/Original_Mess/OriginalMess_SZBF.xlsx')

len = size(DominantMess2);

count = ones(1,1240);

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
writetable(MissmatchedMess2,strcat('SZBF_MissmatchedMess',num2str(libnumber),'.csv'))
end

dd = reshape(count,31,[]);

figure('Renderer', 'painters', 'Position', [10 10 1200 1500])


%h = heatmap(dd,'ColorLimits',[0 1],'GridVisible','off');
h = heatmap(dd','ColorLimits',[0 1]);
h.Title = 'oligo comparison ';
saveas(gcf,strcat('lib',num2str(libnumber),'_SZBF_oligocomparison.jpg'))
writematrix(dd,strcat('lib',num2str(libnumber),'_SZBF_oligocomparison.csv'))

