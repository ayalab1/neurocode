animals = {'Unimplanted\m01','Unimplanted\m02','Unimplanted\m03','Unimplanted\m04','V1Jean','AO52'};
add = [0 0 0 0 17 8];
nDays = 12;
fontsize = 15;

[n_correct,n_right,n_trials,n_skipStart] = deal(nan(length(animals),nDays));
skipStartCell = cell(length(animals),nDays);
for i=1:length(animals)
    animal = animals{i};
    for day = 1:nDays
        folder = fullfile('N:\V1test\',animal,['day' num2str(day + add(i))]);
        try
            [correct,right,skipStart] = PlotYmazePerformance(folder);
            n_correct(i,day) = sum(correct);
            n_right(i,day) = sum(right);
            n_skipStart(i,day) = sum(skipStart>0);
            n_trials(i,day) = length(correct);
            skipStartCell{i,day} = skipStart;
        end
    end
end

%%

clf
p = (n_correct./n_trials); 
p = p*100;

z = nan(1,nDays);
for i=1:nDays,z(i) = zBinomialComparison(sum(n_correct(:,i)),sum(n_trials(:,i)),0.5); end


clf

handles{1} = subplot(2,3,1);
PlotColorMap(p);
set(gca,'ytick',1:length(animals),'yticklabel',out2(@fileparts,animals));
xlabel('training day');
ylabel('animal');
clabel('% correct');
clim([-1 1]*20+50);
set(gca,'FontSize',fontsize);
posAxes{1} = get(gca,'Position');

handles{2} = subplot(2,3,4);
semplot(1:nDays,p);
ylim([0 100]);
ylabel('% correct (mean +/- sem per animal)');
xlabel('training day');
PlotHVLines(50,'h','k--');
sig = FindInterval(z>p2z(0.05));
if ~isempty(sig), sig(:,1) = sig(:,1)-0.5; sig(:,2) = sig(:,2)+0.5; PlotIntervals(sig,'color','r'); end
xlim([1 nDays]);
set(gca,'FontSize',fontsize);
posAxes{2} = get(gca,'Position');

p = (n_right./n_trials); 
p = p*100;
% p(mean(p,2)<50,:) = 100 - p(mean(p,2)<50,:);

z = nan(1,nDays);
for i=1:nDays,z(i) = zBinomialComparison(sum(p(:,i).*n_trials(:,i)/100),sum(n_trials(:,i)),0.5); end


handles{3} = subplot(2,3,2);
PlotColorMap(p);
set(gca,'ytick',1:length(animals),'yticklabel',out2(@fileparts,animals));
xlabel('training day');
ylabel('animal');
clabel('% right (and not left)');
clim([-1 1]*20+50);
set(gca,'FontSize',fontsize);
posAxes{3} = get(gca,'Position');

handles{4} = subplot(2,3,5);
semplot(1:nDays,p);
ylim([0 100]);
ylabel({'% "right" (and not left) trials', '(mean +/- sem per animal)'});
xlabel('training day');
PlotHVLines(50,'h','k--');
sig = FindInterval(z>p2z(0.05));
if ~isempty(sig), sig(:,1) = sig(:,1)-0.5; sig(:,2) = sig(:,2)+0.5; PlotIntervals(sig,'color','r'); end
xlim([1 nDays]);
set(gca,'FontSize',fontsize);
posAxes{4} = get(gca,'Position');


% skip start

n_skipStart = cellfun(@(x) sum(x),skipStartCell);

p = (n_skipStart./n_trials); 
p = p;
% p(mean(p,2)<50,:) = 100 - p(mean(p,2)<50,:);

z = nan(1,nDays);
for i=1:nDays,z(i) = zBinomialComparison(sum(p(:,i).*n_trials(:,i)/100),sum(n_trials(:,i)),0.5); end


handles{5} = subplot(2,3,3);
PlotColorMap(p);
set(gca,'ytick',1:length(animals),'yticklabel',out2(@fileparts,animals));
xlabel('training day');
ylabel('animal');
clabel('Mean # of mouseport skipped');
% clim([-1 1]*20+50);
set(gca,'FontSize',fontsize);
posAxes{5} = get(gca,'Position');

handles{6} = subplot(2,3,6);
semplot(1:nDays,p);
% ylim([0 100]);
ylabel('Mean # of mouseport skipped');
xlabel('training day');
% PlotHVLines(50,'h','k--');
sig = FindInterval(z>p2z(0.05));
if ~isempty(sig), sig(:,1) = sig(:,1)-0.5; sig(:,2) = sig(:,2)+0.5; PlotIntervals(sig,'color','r'); end
xlim([1 nDays]);
set(gca,'FontSize',fontsize);
posAxes{6} = get(gca,'Position');

drawnow

posAxes{2}(3) = posAxes{1}(3); posAxes{4}(3) = posAxes{3}(3); posAxes{6}(3) = posAxes{5}(3);
set(handles{2},'position',posAxes{2});
set(handles{4},'position',posAxes{4});
set(handles{6},'position',posAxes{6});





