function [output]=RNAtoDNA(inputseq)

output = inputseq;

for i = 1:length(inputseq)
    if (output(i) == 'U')
		output(i) = 'T';
	elseif (output(i) == 'u')
		output(i) = 't';
	end
end