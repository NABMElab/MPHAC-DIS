function eff = equaleff(eff1,eff2,cycle)
t1 = 1;
t2 = 1;
for i = 1:cycle
    t1_temp = t1 + t2*eff2;
    t2_temp = t2 + t1*eff1;
    t1 = t1_temp;
    t2 = t2_temp;
end
 t = (t1+t2)/2;
eff = t^(1/cycle);
eff = eff - 1;