close all

figure('Renderer', 'painters', 'Position', [10 10 900 600])

DM = Dominant_count;
DMr = log10(Dominant_count);
plot(DMr,'r.','markersize',20);
xlabel('Different oligos contain 12 English characters');
ylabel('log10(Reads)');
title('DDJ reads uniformity')
set(gca,'fontsize',20);
ylim([1 3]);


figure('Renderer', 'painters', 'Position', [10 10 900 600])

DMrS = sort(DMr,2,'descend');
plot(DMrS,'r.','markersize',20);
xlabel('Different oligos contain 6 English characters');
ylabel('log10(Reads)');
title('DDJ reads uniformity by sort')
set(gca,'fontsize',20);
ylim([1 3]);
hold on
yline(2.250,'b-','linewidth',3);% 80% coverage value
yline(2.373,'k-','linewidth',3);% mean value
legend('reads','80% coverage depth','mean coverage depth')