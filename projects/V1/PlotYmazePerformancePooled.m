animals = {'Unimplanted\m01','Unimplanted\m02','Unimplanted\m03','Unimplanted\m04','V1Jean','AO52'}; % Creates a cell with all the animals.
add = [0 0 0 0 17 8];% Animals sessions have different numbers for the same day, so adding 17 and 8 to V1Jean and AO56 put them in the same number.
fontsize = 15;% Just initiating some parameters.
listOfDays = 20:24; % Specify the days of interest

[n_correct,n_right,n_trials,n_skipStart] = deal(nan(length(animals),length(listOfDays)));% Creates four matrix made out of NaN (not a number). "length(20:nDays)" if I only want to represent the last 4 days, for example.
skipStartCell = cell(length(animals),length(listOfDays));% Creates a cell
for i=1:length(animals) % This looks for all the animals going through all of them
    animal = animals{i};
    for dayNumber = 1:length(listOfDays) % and through all the days/sessions
        day = listOfDays(dayNumber);
        folder = fullfile('N:\V1test\',animal,['day' num2str(day + add(i))]);  % creates a name for each session for each animal
        try % takes the values to the matrix that were created earlier
            [correct,right,skipStart] = PlotYmazePerformance(folder); % "PlotYmazePerformance" is another function that Raly has prepared and it loads the values from the task files.
            n_correct(i,dayNumber) = sum(correct);
            n_right(i,dayNumber) = sum(right);
            n_skipStart(i,dayNumber) = sum(skipStart>0);
            n_trials(i,dayNumber) = length(correct);
            skipStartCell{i,dayNumber} = skipStart;
        end
    end
end

%%

clf
p = (n_correct./n_trials); 
p = p*100; %proportion of correct trials over the total number

z = nan(1,nDays);
for dayNumber = 1:length(listOfDays),z(dayNumber) = zBinomialComparison(sum(n_correct(:,dayNumber)),sum(n_trials(:,dayNumber)),0.5); end


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

z = nan(1,nDays);
for i=1:length(listOfDays),z(i) = zBinomialComparison(sum(p(:,i).*n_trials(:,i)/100),sum(n_trials(:,i)),0.5); end


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

z = nan(1,nDays);
for i=1:length(listOfDays),z(i) = zBinomialComparison(sum(p(:,i).*n_trials(:,i)/100),sum(n_trials(:,i)),0.5); end


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





