function [amplicon] = simPCR(target,primer,cycle,temp,dG,kn,Lt,Lp,Cenz)

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-30);

Template1 = target;
Template2 = target;
Primer1 = primer;
Primer2 = primer;
kf1 = 1e6;
kf2 = 1e6;
kf3 = 1e6;


for i = 1:cycle
    
    [timepoints, concentrations] = ode23s(@preanneal, ...
        linspace(0,30,31), [Template1, Template2, 0, 1e6], options);
    Template1 = concentrations(31,1);
    Template2 = concentrations(31,2);
    dsDNA = concentrations(31,3);
    
    [timepoints, concentrations2] = ode23s(@anneal, ...
        linspace(0,30,31), [Template1, Template2, Primer1, Primer2, dsDNA, 0, 0, dG, temp, kf1, kf2, kf3], options);
    hsDNA1 = concentrations2(31,6);
    hsDNA2 = concentrations2(31,7);
    dsDNA = concentrations2(31,5);
    
    
    [timepoints, concentrations3] = ode23s(@extension, ...
        linspace(0,30,31), [dsDNA, hsDNA1, hsDNA2, kn, Lt, Lp, Cenz], options);
    
    Template1 = concentrations2(31,1)+concentrations3(31,1)+concentrations3(31,2);
    Template2 = concentrations2(31,2)+concentrations3(31,1)+concentrations3(31,3);
    Primer1 = concentrations2(31,3)+concentrations3(31,2);
    Primer2 = concentrations2(31,4)+concentrations3(31,3);
    amplicon(i) = min(Template1,Template2);

    
end 