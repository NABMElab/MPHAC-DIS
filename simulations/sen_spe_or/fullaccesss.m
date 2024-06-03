function needreads = fullaccesss(originalreads, filenumber, efficiency, PCRcycle, needcoverage)

file1 = 2^PCRcycle;
file2 = efficiency^PCRcycle;
file1reads = originalreads*needcoverage;
file2reads = file1reads*(file1/file2)*(filenumber-1);
needreads = file1reads+file2reads;
end

