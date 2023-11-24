% In this tutorial, we will perform Hierarchical Bootstrap.

%% Synthetic data example:

% let's imagine we're comparing the quality of replay in sleep
% following two kinds of tasks:
% Two animals performed task1 as well as task2 and we have 
% a measure of replay quality for each event:

% set up the random number generator so that the results are always the
% same (the "random" numbers will be the same in each execution)
rng(0); 
figure(1);
set(gcf,'position',[1 40 1920 964])% full screen
clf

trueEffect = 0; % Here we set the real effect. 

% This is a synthetic example, in which there will be NO EFFECT of the task
% on replay quality. All the events detected in animal1 have replay scores
% around 0.5, and all the events in animal2 have scores around 0.8, regardless
% of the task. This difference in replay quality may be due to e.g. having
% better recording quality (more cells) in animal2. Here we will produce these
% data:

% Replay quality was around 0.5 for each task in animal1:
animal1task1 = randn(111,1)*0.2 + 0+trueEffect;
animal1task2 = randn(134,1)*0.2 + 0; 

% Replay quality was around 0.8 for each task in animal2:
animal2task1 = randn(567,1)*0.2 + 1+trueEffect; % note that we have a lot more events in sleep following this task: maybe the recording was longer
animal2task2 = randn(200,1)*0.2 + 1; % note that we have very few events in sleep following this task: maybe the animal didn't sleep much

% Plot data:
subplot(2,2,1);
[h,ht] = hist(animal1task1,linspace(-2,2 + trueEffect,100));
[h(2,:),ht] = hist(animal1task2,ht);
[h(3,:),ht] = hist(animal2task1,ht);
[h(4,:),ht] = hist(animal2task2,ht);
plot(ht,Smooth(h(1:2,:)',[2 0]),'linewidth',2);
hold all
set(gca,'ColorOrderIndex',1); % set color back to 1
plot(ht,Smooth(h(3:4,:)',[2 0]),'--','linewidth',2);
set(gca,'ytick',[]); % yaxis is meaningless in a probability distribution
xlabel('replay quality');
legend('task1 (animal1)','task2 (animal1)','task1 (animal2)','task2 (animal2)','box','off','fontsize',15,'location','northwest');
ylim(ylim+[0 diff(ylim)*0.2]);
ylabel('Distribution');
set(gca,'box','off','TickDir','out','fontsize',15);
title('Pooled data distribution');

% If we were to pool these data and perform statistical tests on the pooled data,
% this is what we would get:
subplot(2,2,2);
data = Group([animal1task1;animal2task1],[animal1task2;animal2task2],'sort',false);
anovabox(data(:,1),data(:,2),'alpha',[0 0.05],'boxwidth',0.2); % set the alpha level to zero so that no stars are shown comparing each group to zero. The between-group comparison alpha is 0.05, so if the two groups are different at at least p=0.05, stars will be shown
xlim([0 3]);
pPooled = ranksum([animal1task1;animal2task1],[animal1task2;animal2task2]);
ylabel('replay quality in sleep following task');
set(gca,'xtick',1:2,'xticklabel',{'task1','task2'});
set(gca,'box','off','TickDir','out','fontsize',15);
title(['p(task1 vs task2) = ' num2str(pPooled)]);
title({'Pooled data (not hierarchical):',['p=' num2str(pPooled)]});

% The difference between replay quality in sleep following task1 and task2
% appears to be very significant!! p<10^(-200). Yet remember how in the 
% synthetic data, there is no difference between the two tasks, all these
% differences are because animal2, in which replay quality was better, 
% happened to sleep more after task2 than task1!!

% here we can add the mean replay quality after every task for each animal.
hold on
means = [mean(animal1task1) mean(animal2task1); mean(animal1task2) mean(animal2task2)];
handle = plot(means,'linewidth',2);
legend(handle,{'animal1','animal2'},'box','off','location','southwest');

% We can see that most of the effect is no longer visible.
% Still, in both animals, the mean replay is lower following task2 than task1
% (both lines are moving down slightly).
% Here, with synthetic data, we know that this is due to chance.
% Try changing rng(0) to rng(5), the lines would no longer be going down.
% But with real, nonsynthetic data, we may wonder if this (small) effect
% may be significant. How would we do that, without using a single point for
% each session, which would dramatically reduce our statistical power?

% group data in long format: keep the values at the start and append the animalID
dataTask1 = Group(animal1task1,animal2task1); 
% each row of "dataTask1" contains replay scores for task1 and animalIDs.
% For example, dataTask1(1,:) is [0.1075 1] meanining that a replay score of 
% 0.1075 was observed in animal1.
dataTask2 = Group(animal1task2,animal2task2);
% each row of "dataTask2" contains replay scores for task1 and animalIDs.
% For example, dataTask2(end,:) is [1.0539 2] meanining that a replay score of 
% 1.0539 was observed in animal2.
data = Group(dataTask1,dataTask2);
% each row of "dataTask" contains replay scores, animalIDs, and taskID
% For example, dataTask(1,:) is [0.1075 1 1] meanining that a replay score of 
% 0.1075 was observed in animal1 for task1.
% dataTask(end,:) is [1.0539 2 2] meanining that a replay score of 
% 1.0539 was observed in animal2 for task2.

% "data" is in the correct format to call "HierarchicalBootstrap". Each row is
% a different observation. The first column contains the values that we are
% testing. In this case, replay quality. 
% The next column contains some grouping variable. In this case, the 
% animal ID from which the data comes from. Grouping variables should be 
% ordered from small (e.g. cellID -> session ID -> animal ID, etc). 
% The final column contains the ID of the groups that we want to test.
% In this case, task1 vs task2.

% Run 1000 iterations for the bootstrap, and for each iteration, take the mean replay value for each group:
[p,bootstrapMeans] = HierarchicalBootstrap(data,'nIterations',1000,'fun',@mean); % feel free to try the bootsrap with a different function, e.g. @median

% The way the bootstrap works is that is takes random samples of the data.
% First, it samples on the highest level provided: here, this is the level
% of the animal. For task1, we have provided data from animal1 and animal2.
% The bootstrap will then take a random sample of these two animals.
% This random sample can be (animal1,animal2) or (animal2,animal1) or (animal1,
% animal1) or (animal2,animal2). With more animals, the number of possible
% combinations will be larger. 
% This tries to estimate, if the animals "out there" are as variable as the 
% animals we recorded, then if we picked a random sample of N=2 animals to
% record, we can expect some variablility. Sometimes we'll record 2 of the
% best possible animals, sometimes we'll record 2 of the animals with worst 
% signal we've seen, and most of the time we record a little bit of everything.
% After the animals (groups higher in the hierarchy) have been sampled this way,
% the values within each animal are sampled in the same way (sampling with
% replacement e.g. exactly 111 values from the 111 provided values; some of 
% the original values will be repeated). 

subplot(2,2,3);
% Plot the bootstrapped means:
lims = [floor(min(data(:,1))) ceil(max(data(:,1)))];
[h,ht] = hist(bootstrapMeans,linspace(lims(1),lims(2),50));
bar(ht,h,'edgecolor','none');
% set(gca,'ytick',[]); % yaxis is meaningless in a probability distribution
xlabel('estimated mean replay quality');
legend('task1','task2','box','off','fontsize',15);
ylim(ylim+[0 diff(ylim)*0.2]);
ylabel('number of iterations');
% ylabel('Distribution');
title('Hierarchical bootstrap');
set(gca,'box','off','TickDir','out','fontsize',15);
% Notice how the bootsrapped means of replay quality in task1 (blue) form 3 distinct
% clusters. On the far left, we see a cluster centered on replay quality of zero.
% These are the bootstrap iterations in which we happened to take animal1
% twice (the mean replay quality in animal1 is 0) for task1. Note that
% around a quarter of the iterations (N=250) produced such a result. 
% On the far right, we see a cluster centered on 1. These are the bootstrap 
% iterations in which we happened to take the values of animal2 for task1 twice 
% (the mean replay quality in animal2 is 1). This happened in around a quarter
% of all bootstrap iterations (N=250). At around 0.84, we have a cluster in between.
% These are the bootsrap iterations in which the values for each animal
% were taken once. Animal1 had 111 recorded events (mean=0) and animal2
% had 567 recorded events (mean=1), so the mean of the pooled events in those
% iterations is 0.84. This cluster contains around (N=500 iterations).
% The bootsrapped means of replay quality in task2 (red) form similar
% clusters. The far left and far right clusters are the identical to the 
% blue clusters, because the data within each animal was the same following 
% each task. Only the middle cluster is different, because when pooling data 
% from animal1 as well as animal2, we take 134 'bad' (mean quality = 0) events 
% from  animal1 and 200 'good' (quality=1) events from animal2, respectively, 
% yielding a pooled mean value of around 0.6.
hold all
% One arrow from left to right with text on left side
actualPosition = get(gca,'Position'); actualPosition = [actualPosition(1) actualPosition(3)+actualPosition(1) actualPosition(2) actualPosition(2)+actualPosition(4)];
annotation('textarrow',[1 1]*interp1(xlim',actualPosition(1:2)',0),[0.4 0.3],'String',' animal1 twice ','FontSize',13,'Linewidth',2); % put an arrow at "0"
annotation('textarrow',[1 1]*interp1(xlim',actualPosition(1:2)',1),[0.35 0.25],'String','         animal2 twice ','FontSize',13,'Linewidth',2) % put an arrow at "1"

% The bootstrapped values can be used to produce a p-value:
% in each iteration, we can see the difference between mean replay
% quality in each group (task1 and task2). We plot this difference here:
subplot(2,2,4);
% Plot results
d = -diff(bootstrapMeans,[],2);
lims = [-1 1]*max(abs([min(d) max(d)])); % make symmetrical 
hist(d,linspace(lims(1),lims(2),100));
ylim(ylim+[0 diff(ylim)*0.2]);
PlotHVLines(0,'v','k--','linewidth',2);
xlabel('bootstrapped difference between task1 and task2');
ylabel('number of iterations');
set(gca,'box','off','TickDir','out','fontsize',15);
title(['Hierarchical bootstrap: p=' num2str(p)]);
% Note that here, the distribution is bimodal. The difference is smaller
% (very close to zero, its exact position changes depending on the random
% number generator) when we take the same animal twice, because within the 
% same animal, there is no real difference between replay following the
% two tasks. The other part of the bimodal distribution represents the
% iterations in which each animal is taken once: then, task1 appears to
% result in higher replay quality than task2.
% The p-value is the proportion of the iterations in which the mean
% values of group 1 (replay quality following task one) are lower
% than the mean values of group 2. If the replay quality following
% task1 was better than the replay quality following task2, we expect
% p<0.05.

% This is the end of the tutorial! In this example, the 
% true effect between the two groups was 0 
% This is because we set trueEffect = 0; at the start of the tutorial. 
% Feel free to go back and change the effect size and see how the bootstrap 
% performs when there is an effect, albeit small:
% for example, change it to "trueEffect = 0.05;" and execute the script again.

%% Real data example (data subject to change):
% As an example, we will be answering the question whether
% CA2 barrage pyramidal cells are more likely than chance 
% to connect to a particular subgroup of CA2 interneurons.

% This is the list of sessions we will be using
for i=1
sessionList = {'Y:\Data\SMproject\AO10\day20'
'Y:\Data\SMproject\AO10\day27'
'Y:\Data\SMproject\AO11\day15'
'Y:\Data\SMproject\AO12\day8'
'Y:\Data\SMproject\AO12\day13'
'Y:\Data\SMproject\AO12\day18'
'Y:\Data\SMproject\AO13\day14'
'Y:\Data\SMproject\AO13\day15'
'Y:\Data\SMproject\AO16\day12'
'Y:\Data\SMproject\AO17\day7'
'Y:\Data\SMproject\AO20\day15'
'Y:\Data\SMproject\AO22\day12'
'Y:\Data\SMproject\AO24\day10'
'Y:\Data\SMproject\AO25\day10'
'Y:\Data\SMproject\AO26\day10'
'X:\V1test\AO52\day3'
'X:\V1test\AO52\day6'
'Y:\Data\SMproject\AZ1\day13'
'X:\data\SocialBehavior\SM40\day03'
'X:\data\Barrage\LAK1\day12'
'X:\data\Barrage\NN2\day10'
'X:\data\Barrage\NN2\day10_3'
'X:\data\Barrage\NN2\day15'
'X:\data\Barrage\NN2\day24'
'X:\data\Barrage\NN2\day27'
'X:\data\Barrage\NN2\day28'
'X:\data\Barrage\NN2\day29'
'X:\data\Barrage\NN2\day30'
'X:\data\Barrage\NN2\day31'
'X:\data\Barrage\NN2\day34'
'X:\data\Barrage\NN2\day35'
'X:\data\Barrage\NN2\day38'
'X:\data\Barrage\NN3\day15'
'X:\data\Barrage\NN3\day16'
'X:\data\Barrage\NN3\day25'
'X:\data\Barrage\NN2\day12'
'X:\data\Barrage\NN2\day13'
'X:\data\Barrage\NN2\day14'
'X:\data\Barrage\NN2\day17'
'X:\data\Barrage\NN2\day20'
'X:\data\Barrage\NN2\day21'
'X:\data\Barrage\NN2\day22'
'X:\data\Barrage\NN2\day23'
'X:\data\Barrage\NN2\day81'
'X:\data\Barrage\NN2\day83'
'X:\data\Barrage\NN2\day84'};
end

for ses=1:size(sessionList,1)
    basepath = sessionList{ses};
    [~,animal] = fileparts(fileparts(basepath));
    animals{ses,1} = animal;
end
[~,~,animalID] = unique(animals); % for each session, note the ID of the animal from 1 to N

% For every session, load all the information we need:
dataCell = cell(size(sessionList,1),1);
for ses=1:size(sessionList,1)
    basepath = sessionList{ses};

    disp([datestr(clock) ': Loading data for sesssion ' basepath]);

    cell_metrics = getStruct(basepath,'cell_metrics');

    % First we retrive the list of all CA2 pyramical cells:
    pyr = cellfun(@(x) contains(x,'Pyramidal'), cell_metrics.putativeCellType)';
    ca2 = cellfun(@(x) contains(lower(x),'ca2'), cell_metrics.brainRegion)';
    ca2pyrCells = find(ca2 & pyr);

    % Now, find the CA2 barrage cells in particular:
    load(fullfile(basepath,'Barrage_Files',[basenameFromBasepath(basepath) '.useSpk.UIDkeep.mat']),'UIDkeep');
    ca2barrageCells = UIDkeep(:);

    % The control cells for this analysis will be CA2 pyramidal cells that were not in the CA2 barrage cell list:
    ca2controlCells = ca2pyrCells(~ismember(ca2pyrCells,ca2barrageCells));

    numbers(ses,1:2) = [length(ca2barrageCells) length(ca2controlCells)];
    % Locate the interneurons of interest:
    targetCells = zeros(0,1);
    if isfield(cell_metrics.tags,'Pb')
        targetCells = cell_metrics.tags.Pb(:);
    end
    % make sure they are interneurons:
    targetCells(ismember(targetCells,find(pyr))) = [];

    numbers(ses,1:3) = [length(ca2barrageCells) length(ca2controlCells) length(targetCells)];

    % Make a list of all the possible connections between CA2 pyramidal cells and the interneurons of interest
    [x,y] = ndgrid(ca2pyrCells,targetCells);
    allPossiblePairs = [x(:) y(:)];
    % make sure none of the pairs are between a cell and itself:
    allPossiblePairs(diff(allPossiblePairs,[],2)==0,:) = [];

    % Now we will check, for each of those possible monosynaptic pairs, which ones are actual monosynaptic connections:

    % load the monosynaptic data:
    mono_res = getStruct(basepath,'mono_res');
    monosynaptic = mono_res.sig_con_excitatory;
    monosynaptic(~ismember(monosynaptic(:,2),targetCells),:) = [];
    if isempty(monosynaptic)
        dataCell{ses,1} = zeros(0,4);
    end

    % check for every pair if it is present as a monosynaptic pair
    confirmed = ismember(allPossiblePairs(:,1:2),monosynaptic,'rows');

    % for every pair, check if it's a barrage pair (1) or a control pair (2):
    test = ismember(allPossiblePairs(:,1),ca2barrageCells(:,1));
    data = double(confirmed);
    data(:,2) = ses;
    data(:,3) = animalID(ses);
    data(:,4) = test+1; % so "test" rows will be marked "1" and "control" rows will be marked "2"

    % save these data across sessions:
    dataCell{ses,1} = data;
end


% this gives us the data will all the hierarchical information we need:
% In the first column, we have the actual measure we are testing for. 
% (In this case, the probability for a monosynaptic connection.)
% In the following columns, we have the hierarchical group data 
% (In this case, the session ID in column 3, and the animal ID in column 4).
% The assumtion is that maybe the data within a group are likely to 
% be similiar due to sampling differences (e.g. due to electrode placement,
% animal behavior, etc);
% In the last column, we have the id of the conditions we are comparing.
% (In this case, ca2 barrage cells vs other, control ca2 pyramidal cells).

data = cell2mat(dataCell);
[p,bootstrapMeans] = HierarchicalBootstrap(data,'nIterations',10000);

%%


figure(2); 
set(gcf,'position',[1 40 1920 964])% full screen
clf
% Plot data:
% If we were to pool these data and perform statistical tests on the pooled data,
% this is what we would get:
subplot(2,2,1);
anovabar(data(:,1),data(:,end),'alpha',[0 0.05],'barwidth',0.2);
pPooled = ranksum(data(data(:,end)==1),data(data(:,end)==2),'tail','right');
ylabel('mean connection probability');
title({'Pooled data (not hierarchical):',['p=' num2str(pPooled)]});
set(gca,'xtick',1:2,'xticklabel',{'CA2 pyr barrage cells','other CA2 pyr cells'});
set(gca,'box','off','TickDir','out','fontsize',15);
% If we were to take the mean for each session and perform statistics on this
% per-session data, this is what we would get:
means = Accumulate(data(:,[2 4]),data(:,1),'mode','mean');
subplot(2,2,2);
anovabox(nani(means),[],'alpha',[0 0.05],'boxwidth',0.2);% set the alpha level to zero so that no stars are shown comparing each group to zero. The between-group comparison alpha is 0.05, so if the two groups are different at at least p=0.05, stars will be shown
xlim([0 3]);
hold all
plot(means','linewidth',2);
means = means(~any(isnan(means),2),:);
pPooled = signrank(means(:,1),means(:,2));
ylabel('mean connection probability');
set(gca,'xtick',1:2,'xticklabel',{'CA2 pyr barrage cells','other CA2 pyr cells'});
set(gca,'box','off','TickDir','out','fontsize',15);
title(['p(task1 vs task2) = ' num2str(pPooled)]);
title({'Per-session data (not hierarchical):',['p=' num2str(pPooled)]});

subplot(2,2,3);
% Plot results
[h,ht] = Dist(1000,means);
plot(ht,Smooth(h,[10 0]),'linewidth',2);
set(gca,'ytick',[]); % yaxis is meaningless in a probability distribution
xlabel('estimated mean connection probability');
set(gca,'xtick',1:2,'xticklabel',{'CA2 pyr barrage cells','other CA2 pyr cells'});
ylim(ylim+[0 diff(ylim)*0.2]);
ylabel('Distribution');
title(['p=' num2str(p)]);
set(gca,'box','off','TickDir','out','fontsize',15);

subplot(2,2,3);
% Plot the bootstrapped means:
lims = [(min(bootstrapMeans(:))) (max(bootstrapMeans(:)))];
[h,ht] = hist(bootstrapMeans,linspace(lims(1),lims(2),100));
xlim(quantile(bootstrapMeans(:),[0.001 0.999]));
bar(ht,h,'edgecolor','none');
% set(gca,'ytick',[]); % yaxis is meaningless in a probability distribution
xlabel('estimated mean connectivity');
legend('CA2 pyr barrage cells','other CA2 pyr cells','box','off','fontsize',15);
ylim(ylim+[0 diff(ylim)*0.2]);
ylabel('number of iterations');
title('Hierarchical bootstrap');
set(gca,'box','off','TickDir','out','fontsize',15);

subplot(2,2,4);
d = -diff(bootstrapMeans,[],2);
lims = [-1 1]*max(abs([quantile(d,[0 1])])); % make symmetrical 
hist(d,linspace(lims(1),lims(2),1000));
ylim(ylim+[0 diff(ylim)*0.2]);
PlotHVLines(0,'v','k--','linewidth',2);
xlabel('bootstrapped difference b/n barrage cells and other cells');
xlim(quantile(d,[0.001 0.999]));
ylabel('number of iterations');
set(gca,'box','off','TickDir','out','fontsize',15);
title(['Hierarchical bootstrap: p=' num2str(p)]);








