

target = 1e-15;
primer = 4e-7;
temp = 60;

Template1 = target;
Template2 = target;
Primer1 = primer;
Primer2 = primer;
kbind = 1e6

%dG = [-12 -10 -8.12 -6.93 -5.77 -4.2]
dG = [-12]
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-30);
parfor i = 1:length(dG)
     [timepoints, concentrations{i}] = ode23s(@anneal, ...
        linspace(0,30,31), [Template1, Template2, Primer1, Primer2, 0, 0, 0, dG(i), temp, kbind, kbind, kbind], options);
end

for i = 1:length(dG)
    hold on
    plot(0:1:30,concentrations{i}(:,1))
end
