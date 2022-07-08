
fontsize = 15;% Just initiating some parameters.
% 
% listOfDays = [34 35 36 37 38 39 40 41 42 43 44 45 46 47 48]; % Specify the days of interest
% nDays = length(dir('N:\V1test\Unimplanted\m11'))-2; % the number of folders inside m07, excluding the first two '.' and '..' fields
listOfDays = 15:30;

[n_correct,n_right,n_trials,n_skipStart] = deal(nan(1,length(listOfDays)));% Create four matrix made of NaN (not a number). The dimensions are the number of animals for rows and the number fo days for columns.
skipStartCell = cell(1,length(listOfDays));% Create a cell number of animals (rows) x number of days (columns)
i=1;
for dayNumber = 1:length(listOfDays)
    day = listOfDays(dayNumber);
    if day<10
        folder = ['N:\V1test\V1JeanBaptiste\V1JeanBaptiste_22060',num2str(day)];  % create a name for each session for each animal
    else
        folder = ['N:\V1test\V1JeanBaptiste\V1JeanBaptiste_2206',num2str(day)];  % create a name for each session for each animal
    end

    try % takes the values to the matrix that were created earlier
        [correct,right,skipStart] = PlotYmazePerformance(folder); % "PlotYmazePerformance" is another function that Raly has prepared and it loads the values from the task files.
        %             midpoint = round(length(correct)/2);
        %             correct = correct(midpoint:end); right = right(midpoint:end); skipStart = skipStart(midpoint:end);
        %             n = length(correct); n = floor(n/2);
        %             correct = correct(1:n); right = right(1:n);
        %             skipStart = skipStart(1:n); % only first half of the session.
        %             correct = correct(n+1:end); right = right(n+1:end);
        %             skipStart = skipStart(n+1:end); %second half of the session
        n_correct(i,dayNumber) = sum(correct);
        n_right(i,dayNumber) = sum(right);
        n_skipStart(i,dayNumber) = sum(skipStart>0);
        n_trials(i,dayNumber) = length(correct);
        skipStartCell{i,dayNumber} = skipStart;
        correctCell{i,dayNumber} = correct';
    catch
        disp(['Problem with session ' folder])
    end
end

%%
% figure;
% clf
% maxNtrials = max(cellfun(@length,correctCell));
% for i=1:size(correctCell,2), for j=1:size(correctCell,1), correctCell{j,i} = double(correctCell{j,i}); correctCell{j,i}(end+1:maxNtrials(i)) = nan; end; end
% i=9; this = cell2mat(correctCell(:,i));
% PlotColorMap(nansmooth(this,[0 2]),~isnan(this));
% clim([0 0.8])

%%
figure(1) % write 'figure(2)' if we want to compare plots for the first half of the session versus the second half.
clf
p = (n_correct./n_trials); 
p = p*100; %proportion of correct trials over the total number
pPlot = p; pPlot(isnan(n_trials)) = 50;
 
z = nan(1,length(listOfDays));
for dayNumber = 1:length(listOfDays),z(dayNumber) = zBinomialComparison(sum(n_correct(:,dayNumber)),sum(n_trials(:,dayNumber)),0.5); end


clf

handles{1} = subplot(2,3,1);
imagesc(pPlot); set(gca,'YDir','reverse')
set(gca,'ytick',1:length(animals),'yticklabel',out2(@fileparts,animals));
xlabel('training day');
ylabel('animal');
clabel('% correct');
clim([-1 1]*20+50);
set(gca,'FontSize',fontsize);
posAxes{1} = get(gca,'Position');
ColorMap(gca,[0 0 1],[1 1 1],[1 0 0]);

handles{2} = subplot(2,3,4);
semplot(listOfDays,p);
ylim([0 100]);
ylabel('% correct (mean +/- sem per animal)');
xlabel('training day');
PlotHVLines(50,'h','k--');
sig = FindInterval(z>p2z(0.05));
if ~isempty(sig), sig = listOfDays(sig); sig(:,1) = sig(:,1)-0.5; sig(:,2) = sig(:,2)+0.5; PlotIntervals(sig,'color','r'); end
xlim(listOfDays([1 end]));
set(gca,'FontSize',fontsize);
posAxes{2} = get(gca,'Position');


p = (n_right./n_trials); 
p = p*100;
% p(mean(p,2)<50,:) = 100 - p(mean(p,2)<50,:);

z = nan(1,length(listOfDays));
for i=1:length(listOfDays),z(i) = zBinomialComparison(sum(p(:,i).*n_trials(:,i)/100),sum(n_trials(:,i)),0.5); end


handles{3} = subplot(2,3,2);
imagesc(p); set(gca,'YDir','reverse')
set(gca,'ytick',1:length(animals),'yticklabel',out2(@fileparts,animals));
xlabel('training day');
ylabel('animal');
clabel('% right (and not left)');
clim([-1 1]*20+50);
set(gca,'FontSize',fontsize);
posAxes{3} = get(gca,'Position');
ColorMap(gca,[0 0 1],[1 1 1],[1 0 0]);

handles{4} = subplot(2,3,5);
semplot(listOfDays,p);
ylim([0 100]);
ylabel({'% "right" (and not left) trials', '(mean +/- sem per animal)'});
xlabel('training day');
PlotHVLines(50,'h','k--');
sig = FindInterval(z>p2z(0.05));
if ~isempty(sig), sig = listOfDays(sig); sig(:,1) = sig(:,1)-0.5; sig(:,2) = sig(:,2)+0.5; PlotIntervals(sig,'color','r'); end
xlim(listOfDays([1 end]));
set(gca,'FontSize',fontsize);
posAxes{4} = get(gca,'Position');

% skip start

n_skipStart = cellfun(@(x) sum(x),skipStartCell);

p = (n_skipStart./n_trials); 
p = p;
% p(mean(p,2)<50,:) = 100 - p(mean(p,2)<50,:);

z = nan(1,length(listOfDays));
for i=1:length(listOfDays),z(i) = zBinomialComparison(sum(p(:,i).*n_trials(:,i)/100),sum(n_trials(:,i)),0.5); end


handles{5} = subplot(2,3,3);
imagesc(p); set(gca,'YDir','reverse')
set(gca,'ytick',1:length(animals),'yticklabel',out2(@fileparts,animals));
xlabel('training day');
ylabel('animal');
clabel('Mean # of mouseport skipped');
% clim([-1 1]*20+50);
set(gca,'FontSize',fontsize);
posAxes{5} = get(gca,'Position');
ColorMap(gca,[1 1 1],[0 0 0]);

handles{6} = subplot(2,3,6);
semplot(listOfDays,p);
% ylim([0 100]);
ylabel('Mean # of mouseport skipped');
xlabel('training day');
% PlotHVLines(50,'h','k--');
sig = FindInterval(z>p2z(0.05));
if ~isempty(sig), sig = listOfDays(sig); sig(:,1) = sig(:,1)-0.5; sig(:,2) = sig(:,2)+0.5; PlotIntervals(sig,'color','r'); end
xlim(listOfDays([1 end]));
set(gca,'FontSize',fontsize);
posAxes{6} = get(gca,'Position');

drawnow

posAxes{2}(3) = posAxes{1}(3); posAxes{4}(3) = posAxes{3}(3); posAxes{6}(3) = posAxes{5}(3);
set(handles{2},'position',posAxes{2});
set(handles{4},'position',posAxes{4});
set(handles{6},'position',posAxes{6});

