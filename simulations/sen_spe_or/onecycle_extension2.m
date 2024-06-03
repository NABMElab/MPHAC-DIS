 function efficiency = onecycle_extension2(dG,dG2)
target = 1e-14;
primer = 4e-7;
temp = 60;
kbind = 1e5;
time = 30;
dG_init = dG;
dG_adpt = 0;
length_T = 200;
length_P = 20;
kn = 60;

T1 = target;
T2 = target;
P1 = primer;
P2 = primer;
kf1 = kbind;
kf2 = kbind;
kf3 = kbind;
kf4 = kbind;
kf5 = kn;

LT1 = 0;
LT2 = 0;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-30);
[timepoints, concentrations] = ode23s(@anneal2, ...
        linspace(0,time,time+1), [T1, T2, LT1, LT2, 0, 0, P1, P2, 0, 0, 0, 0, kf1, kf2, kf3, kf4, kf5, dG_init, dG_adpt, temp, length_T, length_P, dG2], options);
% currval = [T1, T2, LT1, LT2, dsDNA1, dsDNA2, P1, P2, hdDNA1, hdDNA2, hdDNA3, hdDNA4, Enzyme, kf1, kf2, kf3, kf4, kf5, dG_init, dG_adpt, temp, length_T, length_P];

efficiency = min((concentrations(30,1) + concentrations(30,5) + concentrations(30,9) +concentrations(30,3) + concentrations(30,6) + concentrations(30,11)), ...
    (concentrations(30,2) + concentrations(30,5) + concentrations(30,10) +concentrations(30,4) + concentrations(30,6) + concentrations(30,12)))/target
