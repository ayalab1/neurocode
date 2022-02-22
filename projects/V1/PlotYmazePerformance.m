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

files = dir(basepath);
for i=1:length(files)
    filenames{i,1} = files(i).name;
end
filename = fullfile(basepath,filenames{find(cellfun(@(x) ~isempty(strfind(x,'task-ymaze-')), filenames),1,'last')});

data = load(filename);
choice = find(abs(data(:,3))==50);
trials = data(choice,:);
trialID = trials(:,4);
nTrials = size(trials,1);
animalWentRight = trials(:,2)<-25; % -20 =mouse went left, -30 = mouse went right. Now animalWentRight=1 means right and animalWentRight=0 means left
trialType = trials(:,5)-2; % trials(:,5) is 2 for left and 3 for right. Now trialType is 1 for right and 0 for left
correct = trials(:,3)>0; % trials(:,3) is 50 for correct and -50 for error
alternation = [false; diff(animalWentRight)~=0];

choice(end+1) = size(data,1);
skipStart = nan(nTrials,1);
for i=1:nTrials
    chunk = data(choice(i):choice(i+1)-1,:);
    photodetectors = chunk(chunk(:,2)<-19 & chunk(:,2)>-31,2);
    photodetectors = photodetectors([true; diff(photodetectors)~=0],:);
    skipStart(i) = length(photodetectors)-1;
end

clf

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
title(basepath);

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


subplot(2,3,4);
counts = [sum(animalWentRight) nTrials; nTrials/2 nTrials];
CountBar(counts,1)
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
PlotHVLines(50,'h','k--');
ylabel('% "right" (and not "left")');
set(gca,'FontSize',fontsize);
title({'Right bias (sig. difference indicates bias)'});

subplot(2,3,5);
counts = [sum(correct) nTrials; nTrials/2 nTrials];
CountBar(counts,1)
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
ylabel('% correct');
set(gca,'FontSize',fontsize);
title({'Performance (sig. difference indicates learning)'});

subplot(2,3,6);
counts = [sum(alternation) nTrials; nTrials/2 nTrials];
CountBar(counts,1)
PlotHVLines(50,'h','k--');
set(gca,'xticklabel',{'data','chance'});
ylim([0 100]);
ylabel('% alternation');
set(gca,'FontSize',fontsize);
title({'Alternation (sig. difference indicates alternation)'});

if ~exist([fullfile(basepath,[sessionID '_behavioral_performance.fig']))
    SaveFig(fullfile(basepath,[sessionID '_behavioral_performance']));
end

% subplot(2,3,6);
% counts = Accumulate([trialType,1-correct]+1); counts(:,end) = sum(counts,2); % [correctLeft totalLeft; correctRight totalRight]
% CountBar(counts,1)
% PlotHVLines(50,'h','k--');
% set(gca,'xticklabel',{'"left" trials','"right" trials'});
% ylim([0 100]);
% ylabel('% probability correct');
% set(gca,'FontSize',fontsize);
% title({'Performance depending on trial type','(sig. difference indicates bias)'});



