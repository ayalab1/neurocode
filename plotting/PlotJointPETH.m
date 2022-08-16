function PlotJointPETH(PETH1,PETH2,smooth,names,durations)

%See JointPETH. This function calls JointPETH and plots the results
% Plase input the names for the 3 variables in the cell "names".
% Example:
% durations = [-1 1];
% [h1,ht] = PETH(ripples.timestamps(:,1),deltaWaves.timestamps(:,1),'nBins',101,'durations',durations);
% [h2,ht] = PETH(barrages.timestamps(:,1),deltaWaves.timestamps(:,1),'nBins',101,'durations',durations);
% smooth = 2; names = {'ripples','barrages','delta waves'}; 
% PlotJointPETH(h1,h2,smooth,names,durations);

clf
if ~exist('smooth','var') || isempty(smooth), smooth = 0; end
if ~exist('durations','var') || isempty(durations), durations = [-1 1]; end
if ~exist('names','var'), names = {'PETH1','PETH2','event'}; end 

[hs{1},hs{2},hs{3}] = JointPETH(PETH1,PETH2,smooth);
ht = linspace(durations(1),durations(2),size(PETH1,2));

subplot(1,3,1);
i=1; PlotColorMap(hs{i},'x',ht,'y',ht); %axis square; c=clim;
hold on; plot(xlim,xlim,'k--','linewidth',2);
c = clim;
axisPosition =  [0.1300  0.300 0.2 0.4];
set(gca,'position',axisPosition);
axisPosition = get(gca,'Position');
xlabel(['(' names{2} ') time from ' names{3} ' (s)']);
set(gca,'FontSize',12,'box','off');
PlotHVLines(0,'v','w--','linewidth',2); PlotHVLines(0,'h','w--','linewidth',2);

ah = axes('Position',[0 0.2 0 0]);
semplot(ht,PETH2,'k',smooth)
set(ah,'box','off');% axis off
ahAxisPosition = axisPosition;
ahAxisPosition(2) = axisPosition(2)+axisPosition(4); % start above the colormap
ahAxisPosition(4) = ahAxisPosition(4)/4; % a quarter of the colormap's height
set(ah,'Position',ahAxisPosition);
PlotHVLines(0,'v','k--','linewidth',2);
title([names{2} ' rate']);
set(gca,'ytick',[]);
set(gca,'FontSize',12,'box','off');

av = axes('Position',[0.2 0.2 0.2 0]);
semplot(ht,PETH1,'k',smooth)
set(av,'box','off');% axis off
avAxisPosition = axisPosition;
avAxisPosition(3) = avAxisPosition(3)/4; % a quarter of the colormap's height
avAxisPosition(1) = axisPosition(1)-avAxisPosition(3); % start on the left of colormap
set(av,'Position',avAxisPosition);
PlotHVLines(0,'v','k--','linewidth',2);
title([names{1} ' rate']);
view([-90 90])
set(gca,'ytick',[]);
set(gca,'FontSize',12,'box','off');

subplot(1,3,2);
i=2; PlotColorMap(hs{i},'x',ht,'y',ht); %axis square; c=clim;
thisPosition = axisPosition; thisPosition(1) = 0.41;
set(gca,'Position',thisPosition);
clim(c);
title(['pooled : expected']);
hold on; plot(xlim,xlim,'k--','linewidth',2);
xlabel(['(' names{2} ') time from ' names{3} ' (s)']); ylabel(['(' names{1} ') time from ' names{3} ' (s)']);
set(gca,'FontSize',12,'box','off');
PlotHVLines(0,'v','w--','linewidth',2); PlotHVLines(0,'h','w--','linewidth',2);

subplot(1,3,3);
thisPosition = axisPosition; thisPosition(1) = 0.69;
set(gca,'Position',thisPosition);
i=3; PlotColorMap(hs{i},'x',ht,'y',ht); %axis square; c=clim;
clim(max(abs(clim))*[-1 1]);
hold on; plot(xlim,xlim,'k--','linewidth',2);
xlabel(['(' names{2} ') time from ' names{3} ' (s)']); ylabel(['(' names{1} ') time from ' names{3} ' (s)']);
set(gca,'FontSize',12,'box','off');
title('observed - expected');
PlotHVLines(0,'v','w--','linewidth',2); PlotHVLines(0,'h','w--','linewidth',2);

ah2 = axes('Position',[0.8 0.3 0 0]);
corrected = CircularShift(hs{i},ceil(size(hs{i},2)/2)-(1:size(hs{i},2)))';
semplot(ht,corrected);
set(av,'box','off');% axis off
ahAxisPosition = thisPosition;
ahAxisPosition(2) = axisPosition(2)+axisPosition(4)+0.05; % start above the colormap
ahAxisPosition(4) = ahAxisPosition(4)/4; % a quarter of the colormap's height
set(ah2,'Position',ahAxisPosition);
PlotHVLines(0,'v','k--','linewidth',2);
title(['corrected ' names{1} ' response to ' names{2}]);
set(gca,'ytick',[]);
set(gca,'FontSize',12,'box','off');
