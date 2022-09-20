function barSumFigs(HSEmetrics)
close all
CA1ses = find(HSEmetrics.Notes == "CA1");
CA2ses = find(HSEmetrics.Notes == "CA2");
CA3ses = find(HSEmetrics.Notes == "CA3");
PYRses = find(HSEmetrics.Notes == "All pyr");

% CA1ses = CA1ses(start:stop);
% CA2ses = CA2ses(start:stop);
% CA3ses = CA3ses(start:stop);
% PYRses = PYRses(start:stop);

avgISI_tot(1) = mean(HSEmetrics.avgISI(CA1ses))*1000;
avgISI_tot(2) = mean(HSEmetrics.avgISI(CA2ses))*1000;
avgISI_tot(3) = mean(HSEmetrics.avgISI(CA3ses))*1000;
avgISI_tot(4) = mean(HSEmetrics.avgISI(PYRses))*1000;
avgFR_tot(1) = mean(HSEmetrics.avgFR(CA1ses));
avgFR_tot(2) = mean(HSEmetrics.avgFR(CA2ses));
avgFR_tot(3) = mean(HSEmetrics.avgFR(CA3ses));
avgFR_tot(4) = mean(HSEmetrics.avgFR(PYRses));
stdISI_tot(1) = std(HSEmetrics.avgISI(CA1ses))*1000;
stdISI_tot(2) = std(HSEmetrics.avgISI(CA2ses))*1000;
stdISI_tot(3) = std(HSEmetrics.avgISI(CA3ses))*1000;
stdISI_tot(4) = std(HSEmetrics.avgISI(PYRses))*1000;
stdFR_tot(1) = std(HSEmetrics.avgFR(CA1ses));
stdFR_tot(2) = std(HSEmetrics.avgFR(CA2ses));
stdFR_tot(3) = std(HSEmetrics.avgFR(CA3ses));
stdFR_tot(4) = std(HSEmetrics.avgFR(PYRses));

figure(1);
hold on
isi = plot(1:4,avgISI_tot, '.');
set(gca,'xtick',[1:4],'xticklabel',{'CA1';'CA2';'CA3';'Pyr'},'xlim',[0 5])
ylabel('ISI (ms)');
title('Average ISI');
er = errorbar(1:4,avgISI_tot,-1*stdISI_tot,stdISI_tot);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off

figure(2);
hold on
plot(1:4,avgFR_tot, '.');
set(gca,'xtick',[1:4],'xticklabel',{'CA1';'CA2';'CA3';'Pyr'},'xlim',[0 5])
ylabel('FR (#/s)');
title('Average FR');
er = errorbar(1:4,avgFR_tot,-1*stdFR_tot,stdFR_tot);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
end