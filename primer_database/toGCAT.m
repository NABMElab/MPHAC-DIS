function [output]=toGCAT(inputseq)
for i = 1:length(inputseq)
    if (floor(inputseq(i)) == 1)
        output(i) = 'a';
	elseif (floor(inputseq(i)) == 2)
		output(i) = 't';
    elseif (floor(inputseq(i)) == 3)
        output(i) = 'c';
    elseif (floor(inputseq(i)) == 4)
        output(i) = 'g';
	elseif (inputseq(i) == 0.1)
		output(i) = ' ';
	else
		output(i) = '?';
		fprintf('Unexpected nucleotide!\n');
	end
	
    if (inputseq(i) ~= floor(inputseq(i)))
		output(i) = upper(output(i));
	end
end