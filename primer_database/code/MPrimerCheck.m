function [badness] = MPrimerCheck(primer, varargin)

%[badness] = MPrimerCheck(primer, varargin)
% returns the badness (primer dimer likelihood) matrix
%   Badness(p1, p2) = Sum(Loss(p1_sub, p2_sub))
%      If subsequences p1_sub and p2_sub are reverse complementary,
%   Loss(p1_sub, p2_sub) = b_gc^(numGC) * b_at^(numAT) / (LA+dist1) / (LA+dist2)
%
% varargin can include:
%  '-LA' (LA)
%     Length attenuation parameter; badness scaled by (LA + dist), where dist
%     is the number of bases to the 3' end.  Default LA = 1.
%  '-BAT' (b_at)
%     Exponential base for AT pairs; default = 2.
%  '-BGC' (b_gc)
%     Exponential base for GC pairs; default = 4.

%addpath('Nablab_Ecalc');

counter = 1;
LA = 1; %Length attenuation % taq=1; phusion=3
hashmod = 402653189; %DO NOT CHANGE THIS NUMBER!
b_at = 2; b_gc = 4;

while (length(varargin) >= counter)
	if (strcmp(varargin{counter}, '-LA'))
		LA = varargin{counter+1};
		counter = counter + 2;
    elseif (strcmp(varargin{counter}, '-BAT'))
        b_at = varargin{counter+1};
		counter = counter + 2;
    elseif (strcmp(varargin{counter}, '-BGC'))
        b_gc = varargin{counter+1};
		counter = counter + 2;
    else        
		error('Unexpected input token %s found!\n', varargin{counter});
	end
end

numprimers = length(primer);

badness = zeros(numprimers, numprimers);

for i = 1:numprimers
    % Set up hash table for primer
    hash = zeros(1, hashmod);
    for m = (length(primer{i})-3):-1:1
        for n = (m+3):length(primer{i})
            subseq = primer{i}(m:n);
            numAT = sum((subseq == 'a')|(subseq == 't'));
            numGC = sum((subseq == 'g')|(subseq == 'c'));
            score = b_at^numAT * b_gc^numGC;
            dist = length(primer{i}) - n;
            
            hashindex = PrimeNumberHash(subseq, hashmod);
            hash(hashindex) = hash(hashindex) + score / (LA+dist);
        end           
    end

    % Calculate badness with every other later primer
    for j = i:numprimers
        for m = (length(primer{j})-3):-1:1
            for n = (m+3):length(primer{j})
                subseq = revcompseq(primer{j}(m:n));
                dist = length(primer{j}) - n;
                
                hashindex = PrimeNumberHash(subseq, hashmod);
                badness(i, j) = badness(i, j) + hash(hashindex) / (LA+dist);
            end
        end        
    end    
end

%Reflect across the diagonal
for i = 1:numprimers
    for j = 1:(i-1)
        badness(i, j) = badness(j, i);
    end
end

end
