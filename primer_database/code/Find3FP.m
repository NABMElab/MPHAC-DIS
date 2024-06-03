function primers = Find3FP(sequence,index,energy,temperature)
%this is function for finding the nested outer, middle & inner primers at a
%specific location. The expected visible region starts from index. All 3
%primers are at upstream of index. 
%example: 
%Exon3Primers =
%Find3FP('ACACGTACGTACGTACGACTACACGTGTCGTGCAGTCGACTGATCG',30,[11.5 12.2],60);


spacer = 2;
outer = FindFP(sequence,index+spacer,energy,temperature);
inner = FindFP(sequence,outer{2,1}(2)+1,energy,temperature);
%mid = FindFP(sequence,inner{2,1}(1)+8,energy,temperature);
for i = 1:3
    primers{i,1} = outer{i,1};
    %primers{i,2} = mid{i,1};
    primers{i,2} = inner{i,1};
end




