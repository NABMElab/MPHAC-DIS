function [curhash] = PrimeNumberHash(text, hashmod)

rep = fromGCAT(lower(text));

curhash = 0;
for i = 1:length(text)
	curhash = curhash * 4 + rep(i);
	curhash = mod(curhash, hashmod);
end

curhash = curhash + 1;


