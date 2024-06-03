clear all
eff1 = 0.2
eff2 = 1
cycle = 30
t1 = 1;
t2 = 1;
for i = 1:cycle
    t1_temp(i) = t1 + t2*eff2;
    t2_temp(i) = t2 + t1*eff1;
    eff1_temp(i) = t1_temp(i)/t1;
    eff2_temp(i) = t2_temp(i)/t2;
    eff3_temp(i,1) = (t1_temp(i)+t2_temp(i))/(t1+t2);
    t1 = t1_temp(i);
    t2 = t2_temp(i);
end

data = [t1_temp;t2_temp; eff1_temp; eff2_temp]
data = transpose(data)



dG1 = [0:0.05:1]
dG2 = 0.2


for i = 1:length(dG1)
    eff4(i,1) = equaleff(dG1(i),dG2,20);
end