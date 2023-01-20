%% Basic pipeline for reactivation
% This pipeline will guide you through three ways to look at reactivation:
% (1) Explained Variance (EV)
% (2) Cell pair correlations (very similar conceptually to EV)
% (3) Detecting task-related components with PCA/ICA and observing the 
%     of those components in sleep epochs
%
% Unlike analyses that require decoding the position of the animal (replay
% and theta sequences), reactivation analysis does not require any positions.
% All you need are: (a) intervals during the task, which is that the patterns
% we are looking for are based on; (b) intervals in pre-task and post-task sleep,
% where we will be looking for reactivation of the behavioral patterns, and,
% of course, (c) spiking data.
% I recommend using theta cycles during running for (a) and activity bursts
% during sws for (b), but other solutions are possible.

%% Load session data:

basepath = 'Z:\Data\Can\OML22\day6'; % we will be using this session for the sake of this tutorial
basename = basenameFromBasepath(basepath);
spikeStructure = importSpikes('basepath',basepath,'CellType',"Pyramidal Cell",'brainRegion','CA1'); % load CA1 pyramidal cells

% We will need a spikes matrix, which all the functions making up this 
% tutorial will expect:
try
    spikes = Group(spikeStructure.times{:},'sort',true);
catch 
    % if there are any problems executing "Group", this is how to produce the [timestamp id] matrix
    spikes = [];
    for i=1:length(spikeStructure.times)
        timestamps = spikeStructure.times{i};
        id = ones(size(timestamps))*i; % the id for that unit
        spikes = [spikes; timestamps id]; % add [timestamps id] for this unit to the matrix
    end
    spikes = sortrows(spikes); % sort spikes accoring to their time
end

%% Detect key intervals

% As stated in the note above, this tutorial uses theta cycles during running 
% as the task period defining the patterns. We load them from here:
behavior = getStruct(basepath,'animal.behavior'); % load behavior structure
% We need running epochs, so make sure the behavior file contains them
if ~isfield(behavior,'run') % Here is some example code to detect running epochs:
    error('No "run" epochs in the behavior file.');
end

% Load NREM sleep epochs:
load(fullfile(basepath,[basename '.SleepState.states.mat']),'SleepState');
NREM = bsxfun(@plus,SleepState.ints.NREMstate,[-1 0]);
% We need to separate between pre-task sleep and post-task sleep.
% Here, I define the task epoch using MergePoints:
MergePoints = getStruct(basepath,'MergePoints');
taskID = find(cellfun(@(x) ~contains(lower(x),'sleep') & ~contains(lower(x),'prob') & ~contains(lower(x),'baseline'),MergePoints.foldernames)); % the task epoch is the epoch that does not mention "sleep" or "prob" or "baseline"
task = MergePoints.timestamps(taskID,:);
% alternatively, you could define it using the running epochs: 
% task = behavior.run, or anything else that makes sense for your project.

preSleep = Restrict(NREM,[0 task(1)]); % any NREM epoch that took place before the start of the first task epoch is pre-task sleep
postSleep = Restrict(NREM,[task(end) Inf]); % any NREM epoch that took place after the end of the last task epoch is pre-task sleep
% Note that this would not work for every project. If you have
% task1 -> sleep -> task2 -> sleep, then the code above would not detect any 
% pre-task sleep. Also, if you task design involves multiple task sessions 
% interleaved by multiple sleep sessions, the simple code above would ignore
% any of the sleep sessions in between the first task and the last task.
% Therefore make sure to define pre-task sleep and post-task sleep in a way
% that makes sense for your project.

% Truncate preSleep and postSleep to one hour:
preSleep = TruncateIntervals(preSleep,3600);
postSleep = TruncateIntervals(postSleep,3600);
% We do this because we know that reactivation (and replay) decreases with 
% sleep. Restricting the analyses to the first hour of NREM sleep considerably
% descreases variability (otherwise it may systematically appear that 
% short sleep sessions have the highest reactivation, and days where you
% happen to record for longer would appear to have poorer reactivation). 

simple = true;

if simple,
    % The simplest possible way to define our intervals would be the following:
    % This is ALTERNATIVE 1:
    taskIntervals = SplitIntervals(behavior.run,'pieceSize',0.02); % 20 ms bins in behavior
    preIntervals = SplitIntervals(preSleep,'pieceSize',0.02); % 20 ms bins in pre-task sleep
    postIntervals = SplitIntervals(postSleep,'pieceSize',0.02); % 20 ms bins in post-task sleep
    % In general, I do recommend restricting the task bins to running periods,
    % because if you include immobility, that also will end up including whatever
    % is being reactivated during awake ripples, which may include remote replay
    % and patterns that are unrelated to the task.
    % Also, I recommend restricting the pre- and post-task intervals to the same
    % brain state, to exclude any effects which might just be due to the
    % awake/sleep composition of the pre- and post- task subsessions.
    % The 20ms timescale is optional, although it is recommended to stay within
    % the 10-100ms range when picking your bin size.
else
    % As an alternative to the simplest solution above, we could define
    % task intervals as theta cycles and pre- and post-task intervals as
    % bursts of activity (during ripples).
    % This is ALTERNATIVE 2:
    try
        thetacycles = getStruct(basepath,'thetacycles'); % Load theta sequences
    catch % if they don't exist, compute them
        error(['No "thetacycles" for .' basepath '. Please run auto_theta_cycles']);
    end
    taskIntervals = Restrict(thetacycles.timestamps,behavior.run);

    % Detect bursts of activity in NREM pre-task and post-task
    % If you have already performed replay analysis, these are already detected:
    try
        replay = getStruct(basepath,'replay'); bursts = replay.intervals;
    catch
        % Otherwise, detect them here for the first time:
        doRestrictToRipples = true; % this is a toggle that decides if you want your bursts to be restricted to detected ripples
        if doRestrictToRipples
            try
                ripples = getStruct(basepath,'ripples');
            catch
                error(['No ripples found for session ' basepath]);
                % If you don't have ripples, you either need to detect them first
                % or set "doRestrictToRipples" to false and use all bursts
            end
        end
        bursts = FindBursts(spikes,'thresholds',[0 2],'smooth',0.01,'durations',[0 1]);
        if doRestrictToRipples
            containsRipple = CountInIntervals(ripples.peaks(:),bursts(:,[1 3]))>0; % find burst that contains ripple in it
            bursts(~containsRipple,:) = []; % discard bursts that do not contain ripples
        end
        % Make sure none of the bursts are theta or gamma bursts associated with running:
        inRunning = IntervalsIntersect(bursts(:,[1 3]),behavior.run);
        bursts(inRunning,:);
    end
    preIntervals = Restrict(bursts(:,[1 3]), preSleep);
    postIntervals = Restrict(bursts(:,[1 3]), postSleep);
end

%% Explained Variance (EV)

% This is the measure introduced by Kudrimoti et al. (1999):
% For a more detailed explanation, read the paper, but the intuitive explanation
% is that the correlation matrix between pairs of cells will change from 
% pre-task sleep to post-task sleep. If this is due to reactivating the task
% experience, we would expect that the correlation matrix constructed during 
% the task will help explain the pre-to-post correlation changes. EV is the 
% percent of the variance from pre-to-post that can be explained using the 
% correlations during the task. As a control, we use the exact same measure
% in the other direction, pretending like pre-task sleep and post-task sleep
% are reversed, yielding the Reverse Explained Variance. If we have more
% task-related activity in post-task sleep than pre-task sleep, we'd expect
% EV>REV. 
[EV, REV] = ExplainedVariance(spikes, preIntervals, taskIntervals, postIntervals);

% Plot the result:
figure(1); clf;
set(gcf,'position',[600 700 900 300]);
bar([EV REV]);
hold all
plot(1:2,[EV REV],'k.-','linewidth',2);
ylim([0 1]);
ylabel('Proportion of variance explained');
set(gca,'xtick',1:2,'xticklabel',{'Explained variance','Reverse explained variance'});
set(gca,'box','off','TickDir','out','fontsize',12);

% Note that this is a single session, and EV was meant to be performed over 
% multiple sessions. Ideally, you will pool all your sessions and perform
% statistical analyses to check if EV is significantly larger than REV.

%% Cell pair correlations

% First, we compute the correlation matrix in behavior:
rTask = CorrInIntervals(spikes, taskIntervals);

% We implement this measure as described by Giri et al (2019):
% The correlation matrices in pre- and post-task sleep are very similar. 
% We can extract the part of the correlation matrix in post-task sleep 
% that cannot be explained by the correlation matrix in pre-task sleep
% (residual correlation):
residuals = CorrResiduals(spikes, preIntervals, postIntervals);
% Note that if the pre- and post- correlations were linearly correlated
% with slope=1 and intercept=0, the the residuals would simply be 
% corrPost - corrPre*1 +0. However, in reality, the best fit might not be
% a one-to-one mapping, perhaps because the firing rates/network excitability
% might systematically increase or decrease all correlations. 

% Now, we see if the residuals of the post-task sleep correlations (what is
% new in post-task sleep after controlling from pre-task sleep) is a good fit
% with the correlations that take place in the task:
[r,p] = nancorr(rTask(:),residuals(:));
% A significant correlation is indication for reactivation.

% Plot the result:
figure(2); clf;
set(gcf,'position',[700 300 500 500]);
plot(rTask(:),residuals(:),'.');
lims = [xlim;ylim]; lims = [min(lims(:,1)) max(lims(:,2))]; xlim(lims); ylim(lims); % make x-axis and y-axis have the same range
axis square
% Fit a line:
ok = ~isnan(rTask(:)) & ~isnan(residuals(:));
a = polyfit(rTask(ok),residuals(ok),1);
hold all
plot(xlim,xlim*a(1) + a(2),'r--','linewidth',2);

ylabel('residual correlation in post-task sleep');
xlabel('task correlation');
title(['r= ' num2str(r) ', p=' num2str(p)]);
set(gca,'box','off','TickDir','out','fontsize',12);


% Here, we have used firing correlations in the task for the x-axis.
% Alternative approaches include using overlap between place fields,
% correlations between firing maps, etc. 

%% Activation strength of task-derived components

% First, detect the activity templates expressed in the task:
rng(0); % Because the fastICA algorithm depends on initialization, its results
% will fluctuate slightly between implementations. This can be annoying, 
% because sometimes it is useful to retrace our steps and re-run the same 
% code, and having slightly different components (e.g. component 4 from the 
% previous run may now be component 6). To make sure we can reproduce our
% results, we set the random number generator with a value (here, 0), which
% will result in the same output every time we execute the following line:
templates = ActivityTemplates(spikes,'bins',taskIntervals,'mode','ica');
% Then, detect the (re)activation of each template/component over the whole session
binSize = 0.02; step = 0.001;
bins = Bins(0,MergePoints.timestamps(end),binSize,step); 
% We don't necessarity need the bins outside of pre-task sleep and post-task sleep
bins = Restrict(bins,[preSleep; postSleep]); % this saves memory when computing the strength
strength = ReactivationStrength(spikes,templates,'bins',bins);
% The first column of "strength" is timestamps and the following columns are the activation strength of each component for that timestamp

% Compute the activation strength around pre-task and post-task ripples:
% If you have not yet detected bursts, choosing the "simple" option above, detect them now:
ripples = getStruct(basepath,'ripples');
preRipples = Restrict(ripples.timestamps(:,1),preSleep);
postRipples = Restrict(ripples.timestamps(:,1),postSleep);
nComponents = size(templates,3);
nBins = 101; durations = [-1 1]*0.5;
[mPre,mPost] = deal(nan(nComponents,nBins)); % pre-allocate matrices
for i=1:nComponents
    % use "mPETH" to compute the mean peri-event time histogram around ripples
    mPre(i,:) = mPETH(strength(:,[1 1+i]),preRipples,'durations',durations,'nBins',nBins);
    mPost(i,:) = mPETH(strength(:,[1 1+i]),postRipples,'durations',durations,'nBins',nBins);
end
x = linspace(durations(1),durations(2),nBins);

%% Plot the result:
figure(3); clf;
set(gcf,'position',[300 200 1500 700]);
matrices = {mPre,mPost,mPost-mPre};
colors = {[0 0 0],[0.5 0 0],[1 0 0]};
titles = {'pre-task sleep','post-task sleep','difference'};
for i=1:3
    subplot(2,4,i); PlotColorMap(matrices{i},'x',x);
    set(gca,'box','off','TickDir','out','fontsize',12);
    ylabel('component ID');
    PlotHVLines(0,'v','w--','linewidth',2);
    title(titles{i});
end
for i=1:3
    subplot(2,4,i+4); 
    semplot(x,matrices{i},colors{i});
    xlabel('time from ripple start (s)');
    ylabel('mean activation strength');
end
clims % set the same color limits on the three panels
EquateScales(4,5,6); % set the same y-axis limits on the bottom three panels

% add a line at 0
for i=1:3
    subplot(2,4,i+4);
    PlotHVLines(0,'v','k--','linewidth',2);
    PlotIntervals([0 0.1]);
    set(gca,'box','off','TickDir','out','fontsize',12);
end

% define the mean response to ripples as the strength between 0 and 100ms after the ripple start:
response = [mean(mPre(:,x>0 & x<0.1),2) mean(mPost(:,x>0 & x<0.1),2)];

subplot(1,4,4);
anovabox(response,[],'alpha',[0 0.05]);
hold on
funColors = true; % if you want to have fun colors:
if funColors
    % make colormap from red (strongest component) to blue (weakest component)
    % (the components are ordered by the variance they explain)
    map = zeros(size(response,1),3);
    map(1,:) = [1 0 0]; % red
    map(end,:) = [0 0 1]; % blue
    for i=1:3, map(:,i) = linspace(map(1,i),map(end,i),size(map,1)); end
    set(gca,'ColorOrder',map,'ColorOrderIndex',1); hold on
end
handle = plot(response','.-','linewidth',1.5);
p = signrank(response(:,1),response(:,2));
title(['difference: p=' num2str(p)]);
set(gca,'xtick',1:2,'xticklabel',{'pre','post'});
ylabel('Response 0-100ms after ripple start');
set(gca,'box','off','TickDir','out','fontsize',12);
names = cell(1,nComponents); 
for i=1:nComponents, names{i} = ''; end
names{1} = 'Strongest component (1)'; names{nComponents} = ['Weakest component (' num2str(nComponents) ')']; 
legend(handle,names);
legend('box','off','location','northoutside');









