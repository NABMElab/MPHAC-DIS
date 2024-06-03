function [output]=revcomp(inputseq)
output = zeros(1,length(inputseq));
for i = 1:length(inputseq)
    if ((floor(inputseq(i)) == 1)||(floor(inputseq(i)) == 3))
        output(i) = inputseq(i)+1;
    elseif ((floor(inputseq(i)) == 2)||(floor(inputseq(i)) == 4))
        output(i) = inputseq(i)-1;
	elseif (inputseq(i) == 0.1)
		output(i) = 0.1;
	else
        fprintf('Unexpected nucleotide!\n');
    end
end
output = fliplr(output);