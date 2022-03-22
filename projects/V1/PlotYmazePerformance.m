function [correct,animalWentRight,skipStart] = PlotYmazePerformance(basepath)

if ~exist('basepath','var')
    basepath = pwd;
end

[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
fontsize = 15;
markersize = 20;
smooth = 2; % for the semplots

z2p = @(z) cdf('norm',-(abs(z)),0,1)*2; % predefine a function converting z-values to p-values

files = dir(basepath);
for i=1:length(files)
    filenames{i,1} = files(i).name;
end
indices = find(cellfun(@(x) ~isempty(strfind(x,'task-ymaze-')), filenames));

data = []; currentTime = 0; currentTrial = 0;
for i=1:length(indices)
    filename = fullfile(basepath,filenames{indices(i)});
    loaded = load(filename);
    loaded(2:end,1) = loaded(2:end,1)+currentTime; 
    nTrialsPerSession(i) = loaded(end,4);
    loaded(:,4) = loaded(:,4)+currentTrial;
    currentTime = loaded(end,1); currentTrial = loaded(end,4);
    data = [data; loaded];
end

choice = find(abs(data(:,3))==50);
trials = data(choice,:);
trialID = trials(:,4);
nTrials = size(trials,1);
animalWentRight = trials(:,2)<-25; % -20 =mouse went left, -30 = mouse went right. Now animalWentRight=1 means right and animalWentRight=0 means left
trialType = trials(:,5)-2; % trials(:,5) is 2 for left and 3 for right. Now trialType is 1 for right and 0 for left
correct = trials(:,3)>0; % trials(:,3) is 50 for correct and -50 for error
alternation = [false; diff(animalWentRight)~=0];
trialDuration = [diff(trials(:,1)); nan];

choice(end+1) = size(data,1);
skipStart = nan(nTrials,1);
for i=1:nTrials
    chunk = data(choice(i):choice(i+1)-1,:);
    photodetectors = chunk(chunk(:,2)<-19 & chunk(:,2)>-31,2);
    photodetectors = photodetectors([true; diff(photodetectors)~=0],:);
    skipStart(i) = length(photodetectors)-1;
end

if exist(fullfile(basepath,[sessionID '_behavioral_performance.fig']),'file')
    return
end

figure(1)
clf
set(gcf,'position',[1 1 1920 1000]); % make fugure big so legends are visible when exported

subplot(2,2,1); cla
colors = get(gca,'ColorOrder');
plot(trialID,animalWentRight,'ko-','markersize',markersize/2,'linewidth',2);
hold all
plot(trialID(correct),animalWentRight(correct),'.','markersize',markersize,'color',colors(1,:));
plot(trialID(~correct),animalWentRight(~correct),'r.','markersize',markersize,'color',colors(2,:));

nans = ones(nTrials,nTrials);
matrix = double(repmat(alternation,1,nTrials));
nans(sub2ind(size(matrix),(1:nTrials)',(1:nTrials)')) = 0;
nans = Smooth(nans,[0 smooth]);
matrix(nans>0.95) = nan;
semplot(1:nTrials,matrix,'y');

nans = ones(nTrials,nTrials);
matrix = double(repmat(animalWentRight,1,nTrials));
nans(sub2ind(size(matrix),(1:nTrials)',(1:nTrials)')) = 0;
nans = Smooth(nans,[0 smooth]);
matrix(nans>0.95) = nan;
semplot(1:nTrials,matrix,'k');
PlotHVLines(0.5,'h','k--','linewidth',2);
repeatSlash = sort([(1:length(basepath)) strfind(basepath,'\')]); % avoid a warning message telling us that the slash was unexpected
title(basepath(repeatSlash));

handle = legend('chosen direction','correct','error','alternation');
set(handle,'location','best','FontSize',fontsize);
set(gca,'ytick',[0 1],'yticklabel',{'left','right'});
xlabel('trial ID');
ylabel('Chosen direction');
set(gca,'FontSize',fontsize);

subplot(2,2,2); cla
colors = get(gca,'ColorOrder');
plot(trialID,correct,'ko-','markersize',markersize/2,'linewidth',2);
hold all
plot(trialID(~animalWentRight),correct(~animalWentRight),'r.','markersize',markersize,'color',colors(4,:));
plot(trialID(animalWentRight),correct(animalWentRight),'.','markersize',markersize,'color',colors(3,:));
PlotHVLines(0.5,'h','k--','linewidth',2);
title(['Number of trials: ' num2str(nTrialsPerSession)]);

nans = ones(nTrials,nTrials);
matrix = double(repmat(correct,1,nTrials));
nans(sub2ind(size(matrix),(1:nTrials)',(1:nTrials)')) = 0;
nans = Smooth(nans,[0 smooth]);
matrix(nans>0.95) = nan;
semplot(1:nTrials,matrix,'k');

handle = legend('correct','left','right');
set(handle,'location','best','FontSize',fontsize);
set(gca,'ytick',[0 1],'yticklabel',{'error','correct'});
xlabel('trial ID');
ylabel('Animal choice');
set(gca,'FontSize',fontsize);

subplot(2,2,3); cla
colors = get(gca,'ColorOrder');
plot(trialID,skipStart,'ko-','markersize',markersize/2,'linewidth',2);
hold all
plot(trialID(~correct),skipStart(~correct),'r.','markersize',markersize,'color',colors(1,:));
plot(trialID(correct),skipStart(correct),'.','markersize',markersize,'color',colors(2,:));
PlotHVLines(0.5,'h','k--','linewidth',2);

nans = ones(nTrials,nTrials);
matrix = double(repmat(skipStart,1,nTrials));
nans(sub2ind(size(matrix),(1:nTrials)',(1:nTrials)')) = 0;
nans = Smooth(nans,[0 smooth]);
matrix(nans>0.95) = nan;
semplot(1:nTrials,matrix,'k');

handle = legend('skipStart','correct','error');
set(handle,'location','best','FontSize',fontsize);
set(gca,'ytick',[0 1],'yticklabel',{'go back to start','skipStart'});
xlabel('trial ID');
ylabel('Animal choice');
set(gca,'FontSize',fontsize);

subplot(2,2,4); cla
colors = get(gca,'ColorOrder');
plot(trialID,trialDuration,'ko-','markersize',markersize/2,'linewidth',2);
hold all
plot(trialID(~correct),trialDuration(~correct),'r.','markersize',markersize,'color',colors(1,:));
plot(trialID(correct),trialDuration(correct),'.','markersize',markersize,'color',colors(2,:));

nans = ones(nTrials,nTrials);
matrix = double(repmat(trialDuration,1,nTrials));
nans(sub2ind(size(matrix),(1:nTrials)',(1:nTrials)')) = 0;
nans = Smooth(nans,[0 smooth]);
matrix(nans>0.95) = nan;
semplot(1:nTrials,matrix,'k');

handle = legend('trial duration','correct','error');
set(handle,'location','best','FontSize',fontsize);
% set(gca,'ytick',[0 1],'yticklabel',{'0','skipStart'});
xlabel('trial ID');
ylabel('Trial duration (s)');
set(gca,'FontSize',fontsize);
drawnow

if true | ~exist(fullfile(basepath,[sessionID '_behavioral_performance.fig']))
    SaveFig(fullfile(basepath,[sessionID '_behavioral_performance']));
end

figure(2)
clf
set(gcf,'position',[1921 1 1920 1000]); % make fugure big so legends are visible when exported

subplot(2,3,1);
counts = [sum(animalWentRight) nTrials 0.5];
BinomialBar(counts)
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
PlotHVLines(50,'h','k--');
ylabel('% "right" (and not "left")');
set(gca,'FontSize',fontsize);
xlabel({'Right bias','(sig. difference indicates bias)'});
z = zBinomialComparison(counts(1),counts(2),counts(3));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');
title(basepath(repeatSlash));

subplot(2,3,2);
counts = [sum(correct) nTrials 0.5];
BinomialBar(counts)
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
ylabel('% correct');
set(gca,'FontSize',fontsize);
title({'Performance','(sig. difference indicates learning)'});
z = zBinomialComparison(counts(1),counts(2),counts(3));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');

subplot(2,6,5);
half = floor(length(correct)/2);
counts = [sum(correct(1:half)) half 0.5];
BinomialBar(counts)
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
ylabel('% correct');
set(gca,'FontSize',fontsize);
title({['First half (trials 1-' num2str(half) ')']});
z = zBinomialComparison(counts(1),counts(2),counts(3));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');

subplot(2,6,6);
counts = [sum(correct(half+1:end)) half 0.5];
BinomialBar(counts)
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
ylabel('% correct');
set(gca,'FontSize',fontsize);
title({['Second half (trials ' num2str(half+1) '-' num2str(nTrials) ')']});
z = zBinomialComparison(counts(1),counts(2),counts(3));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');

subplot(2,4,5);
counts = [sum(alternation) nTrials 0.5];
BinomialBar(counts)
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
ylabel('% alternation');
set(gca,'FontSize',fontsize);
title({'Alternation','(sig. difference indicates alternation)'});
z = zBinomialComparison(counts(1),counts(2),counts(3));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');

subplot(2,4,6);
counts = Accumulate([correct 1-alternation]+1); counts(:,end) = sum(counts,2);
CountBar(counts,1);
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'following error trials','following correct trials'},'XTickLabelRotation',15);
ylim([0 100]);
ylabel('% probability to switch');
set(gca,'FontSize',fontsize);
title({'''switch when wrong'' strategy',' '});
z = zBinomialComparison(counts(1,1),counts(1,2),counts(2,1),counts(2,2));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');

subplot(2,4,7);
try anovabar(trialDuration,correct); catch bar(mean(trialDuration)); if mean(correct)>0.5, set(gca,'xticklabel','correct trials'); else  set(gca,'xticklabel','error trials'); end; end
set(gca,'xticklabel',{'error trials','correct trials'},'XTickLabelRotation',15);
ylabel('mean trial duration (s)');
set(gca,'FontSize',fontsize);
title({'correct when engaged?','(are error trials longer?)'});
z = zBinomialComparison(counts(1,1),counts(1,2),counts(2,1),counts(2,2));
handle = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(handle,['p=' num2str(z2p(z))],'box','off');

subplot(2,4,8); cla
block = ceil((1:nTrials)'/10)*10;
anovabar(trialDuration,block,0.05,false,true); 
ylabel('mean trial duration (s)');
xlabel('trial number (blocks of 10)');
set(gca,'FontSize',fontsize);
title({'Task engagement',' '});

drawnow
if true | ~exist(fullfile(basepath,[sessionID '_behavioral_performance_quantification.fig']))
    SaveFig(fullfile(basepath,[sessionID '_behavioral_performance_quantification']));
end


