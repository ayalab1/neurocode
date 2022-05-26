batch = StartBatch(@BatchCanPhasePrecession,'OMLproject.batch');
X = get(batch,'UserData');

%%
X = get(batch,'UserData');
X(1:11,:) = [];
SpectrogramCell{1} = cell2mat(X(:,1)); SpectrogramCell{2} = cell2mat(X(:,2));
f = X{1,3};

M = [];
interval = [0.5 1.5];
thisF = f(InIntervals(f,interval));
clf
colors = {'r','k'};
clear hmhm
for k=1:2,
    sp = SpectrogramCell{k}(:,InIntervals(f,interval));
    z = zscore(Smooth(sp,[0 0]),[],2);
    [~,m] = max(z,[],2);
    % Remove cells without peaks around theta (+/- 25%)
    bad = ~InIntervals(thisF(m),[0.75 1.25]); 
    z(bad,:) = [];
    id = find(~bad);
    [~,m] = max(z,[],2); [~,o] = sort(m); o = 1:length(o);
    z = z(o,:);
    subplot(2,3,(k-1)*3+1);hold off
    PlotColorMap(z,'x',thisF,'bar','on'); hold on;
    plot(thisF(m(o)),1:length(m),'k.-','markersize',1)
    clabel('power');
    PlotHVLines(1,'v','w--','linewidth',2);
    if k==1;title('stim');elseif k==2, title('z-scored spectra, nonstim');end
    if k==1, ylabel('neuron ID (ordered)'); end
    xlabel('frequency (spikes per theta cycle)');
    clim([-1 1]*2);
    subplot(2,3,(k-1)*3+2); hold off;
    %     semplot(f,z);
    hist(thisF(m),30);
%     hist(thisF(m),10);
    PlotHVLines(1,'v','k--','linewidth',2);
    m = m(:); m(:,2)=k; m(:,3) = id; 
    M = [M;m];
    hmhm(:,k) = double(~bad); hmhm(bad,k) = nan;
    hmhm(~bad,k) = m(:,1);
    if k==1;title('distribution of neurons'' peak firing frequency, stim epochs');else, title('nonstim epochs'); end
    xlabel('Firing frequency (per theta cycle)');
    xlim([0.5 1.5]);

     subplot(2,3,3); if k==1, hold off; end
     semplot(thisF,z,colors{k},1);
     if k==2, PlotHVLines(1,'v','k--','linewidth',2); end

     xlabel('Firing frequency (per theta cycle)');
end
value = thisF(M(:,1))-1;
group = M(:,2);
stats = [out2(@kstest2,value(group==1),value(group==2))];
subplot(2,3,6); cla
anovaplot(value+1,group,'parametric','off','alpha',[0 0.05]);
p = ranksum(value(group==1),value(group==2));
title(['ranksum: stim vs nonstim, p=' num2str(p(1))]);
p = [signrank(value(group==1),[],'tail','right'), signrank(value(group==2),[],'tail','right')];
set(gca,'xtick',1:2,'xticklabel',{['stim>1, p=' num2str(p(1))],['nonstim>1, p=' num2str(p(2))]});
ylabel('Peak firing frequency (per theta cycle)');
subplot(2,3,5);
ylabel(['kstest: stim vs nonstim, p=' num2str(stats(1))]);





