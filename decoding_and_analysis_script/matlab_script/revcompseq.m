function [output]=revcompseq(inputseq)

output = toGCAT(revcomp(fromGCAT(inputseq)));
