function [] = SummarizeBadness(primer, badness, varargin)

%[] = SummarizeBadness(primer, badness)
% returns text report of key contributors to badness
%  '-TOP' (topnum)
%     Number of top predicted primer dimers to display

addpath('Nablab_Ecalc');

counter = 1;
topnum = 3;

while (length(varargin) >= counter)
	if (strcmp(varargin{counter}, '-top'))
		topnum = varargin{counter+1};
		counter = counter + 2;
    else
		error('Unexpected input token %s found!\n', varargin{counter});
	end
end

numprimers = length(primer);

fprintf(['\nNumber of primers: ', num2str(numprimers)]);

for i = 1:numprimers
    plength(i) = length(primer{i});
end

fprintf(['\nPrimer length range: ', num2str(min(plength)), 'nt to ', num2str(max(plength)), 'nt']);

totalbadness = sum(sum(badness));
fprintf(['\nBadness sum = ', num2str(totalbadness)]);

for i = 1:numprimers
    pbadness(i) = sum(badness(i,:));
end

fprintf(['\n\n    Top ', num2str(topnum), ' bad primers:']);

%Worst primers
for i = 1:topnum
    pindex = find(pbadness == max(pbadness));
    fprintf(['\n#', num2str(i), ': Primer ', num2str(pindex)]);
    fprintf(['\n  Sequence: ', primer{pindex}]);
    fprintf(['\n  Badness: ', num2str(pbadness(pindex)), ' (', num2str(100*pbadness(pindex)/totalbadness), ' percent)']);
    pbadness(pindex) = 0;
end

fprintf(['\n\n    Top ', num2str(topnum), ' bad primer pairs:']);

%Worst primer pairs
for i = 1:topnum
    pairindex = find(badness == max(max(badness)));
    p1index = rem(pairindex(1), numprimers);
    if p1index == 0
        p1index = numprimers;
    end
    p2index = floor((pairindex(1)-1)/ numprimers)+1;
    fprintf(['\n#', num2str(i), ': Primers ', num2str(p1index), ' and ', num2str(p2index)]);
    fprintf(['\n  Primer ', num2str(p1index), ': ', primer{p1index}]);
    fprintf(['\n  Primer ', num2str(p2index), ': ', primer{p2index}]);
    if p1index == p2index
        fprintf(['\n  Badness: ', num2str(badness(p1index, p2index)), ' (', num2str(100*badness(p1index,p2index)/totalbadness), ' percent)']);
    else
        fprintf(['\n  Badness: ', num2str(badness(p1index, p2index)), ' (', num2str(2*100*badness(p1index,p2index)/totalbadness), ' percent)']);
    end
    
    badness(p1index, p2index) = 0;
    badness(p2index, p1index) = 0;
    
    % badness
end

fprintf('\n\n');

end
