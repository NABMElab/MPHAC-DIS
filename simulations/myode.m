function [timepoints, outputval] = myode(ratefunction, timepoints, initval)
deltat = 0.01;
outputval = [];
currt = 0;
currval = initval;
currindex = 1;
while (currt < timepoints(end))
rates = feval(ratefunction, currval);
currval = currval + rates * deltat;
currt = currt + deltat;
if (currt > timepoints(currindex))
outputval = [outputval; currval];
currindex = currindex + 1;
end
end
myode(TwoStateHyb,[1:30:1800],[1e7,1e7,0,1e6,0])