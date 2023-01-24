%% Basic pipeline for theta sequences
% These are the basic functions required to do theta sequence analysis
% This is a good test to see if your data is good enough to perform
% replay analysis. If theta sequences don't work, replay is not likely to 
% work.
% Make sure you have clean linearized positions before starting this.

%% Load spikes, positions, and theta cycles

basepath = 'Z:\Data\Can\OML22\day7'; % we will be using this session for the sake of this tutorial
basename = basenameFromBasepath(basepath);
spikeStructure = importSpikes('basepath',basepath,'CellType',"Pyramidal Cell",'brainRegion','CA1'); % load CA1 pyramidal cells
behavior = getStruct(basepath,'animal.behavior'); % load behavior structure

% We will be needing linearized positions:
positions = [behavior.timestamps(:) behavior.position.linearized(:)];
positions(any(isnan(positions),2),:) = [];
positions(:,2) = ZeroToOne(positions(:,2)); 

% We will need a spikes matrix:
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


try
    thetacycles = getStruct(basepath,'thetacycles'); % Load theta sequences
catch % if they don't exist, compute them
    auto_theta_cycles 
    thetacycles = getStruct(basepath,'thetacycles');
end

%%

% We need running epochs, so make sure the behavior file contains them
if ~isfield(behavior,'run') % Here is some example code to detect running epochs:
    % Make sure speed is correctly computed: remove position outliers and recompute speed:
    t = (behavior.timestamps(1):mode(diff(behavior.timestamps(:))):behavior.timestamps(end))';
    ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:));
    interpolated = interp1(behavior.timestamps(ok)',[behavior.timestamps(ok)' behavior.position.x(ok)',behavior.position.y(ok)'],t);
    interpolated(isnan(interpolated(:,1)),:) = [];
    smoothed = behavior.position.x; smoothed(~isnan(smoothed)) = Smooth(smoothed(~isnan(smoothed)),5);
    dx = abs(smoothed-behavior.position.x);
    smoothed = behavior.position.y; smoothed(~isnan(smoothed)) = Smooth(smoothed(~isnan(smoothed)),5);
    dy = abs(smoothed-behavior.position.y);
    d = sqrt(dx.^2 + dy.^2);
    outliers = isoutlier(d,'ThresholdFactor',50);
    behavior.position.x(outliers) = nan; behavior.position.y(outliers) = nan; behavior.position.linearized(outliers) = nan;
    speed = LinearVelocity([behavior.timestamps(:) behavior.position.x(:),behavior.position.y(:)],15); t = speed(:,1);
    % define running as speed above 20% of the peak speed in a trial. Compute trials:
    allPositions = [behavior.timestamps(:) behavior.position.linearized(:)]; allPositions(:,2) = ZeroToOne(allPositions(:,2));
    [up,down] = ZeroCrossings([allPositions(:,1),allPositions(:,2)-0.5]);
    midpoints = allPositions(up|down,1);
    ok = InIntervals(allPositions(FindClosest(allPositions(:,1),midpoints),2),[0.3 0.7]); % true midpoints are close to positions around 0.5 (not 1 to 0 wrap crossings)
    midpoints = midpoints(ok);
    topSpeed = interp1(speed(~any(isnan(speed),2),1),speed(~any(isnan(speed),2),2),midpoints);
    threshold = nanmedian(topSpeed)/5; % the threshold is 20% of the peak speed
    run = t(FindInterval(speed(:,2)>threshold));
    run = ConsolidateIntervals(run,'epsilon',0.5); % merge nearby run epochs
    % Make sure periods without position data are not erroneously included in "running"
    t0 = (behavior.timestamps(1):mode(diff(behavior.timestamps(:))):behavior.timestamps(end))';
    gaps = abs(reshape(behavior.timestamps(FindClosest(behavior.timestamps(:),t0(:))),[],1) - t0(:));
    noData = t0(FindInterval(gaps>0.1)); if size(noData,2)==1, noData=noData'; end
    run = SubtractIntervals(run,noData);
    run = ConsolidateIntervals(run,'epsilon',0.5);
    run(diff(run,[],2)<0.5,:) = []; % remove run epochs lasting for less than 0.5s
    behavior.run = run;
    behavior.speed = speed(:,2)'; % replace old speed with this speed

    % I recommend saving the running epochs in the behavior file:
    save(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');
end

ok = InIntervals(thetacycles.timestamps,behavior.run);
cycles = thetacycles.timestamps(ok,:);
peaks = thetacycles.peaks(ok);
if ~isfield(thetacycles,'amplitude'), amplitude = thetacycles.amplitude(ok); else, amplitude = zeros(0,1); end

%% Perform Bayesian recostruction on the (slow) behavioral timescale

% First, take a look if the reconstructed position is a good fit with the
% actual position. It is good practice to split the dataset in two: one 
% half to train the decoder, and another half to actually decode the 
% position. This way, if you can decode the position well, you know it is 
% because the spikes contain information on the animal position and not 
% overfitting.

nBins = 200; timescale = 0.1;
bins = SplitIntervals(behavior.run,'pieceSize',timescale,'discard',false); % split running epochs into one-second bins
half1 = rem((1:size(bins,1))',2)==0;
half2 = ~half1; 
estimations = nan(nBins,size(bins,1)); actual = nan(size(bins,1),1);

[estimations1,actual1] = ReconstructPosition(positions,spikes,bins(half1,:),'training',bins(half2,:));
estimations(:,half1) = estimations1; actual(half1) = actual1;
[estimations2,actual2] = ReconstructPosition(positions,spikes,bins(half2,:),'training',bins(half1,:));
estimations(:,half2) = estimations2; actual(half2) = actual2;

figure(1); % visually inspect the reconstruction quality
clf;
x = (1:size(estimations,2))'; % the bin number for each estimation
yPosition = linspace(0,1,size(estimations,1));
PlotColorMap(estimations,'x',x,'y',yPosition,'cutoffs',[0 0.01]);
ColorMap(gca,[1 1 1],[0 0 0]); % change color scale to white-to-black
% xlim([0 100]);
hold on
plot(x,actual,'r.-','linewidth',1.5); % plot the true position of the animal
xlim([1500 1700]);
xlabel('Theta cycles')
ylabel('Linearized position');
clabel('Probability (low -> high)');
set(gca,'box','off','TickDir','out','fontsize',12);
% Note how the decoded position tends to fit the true animal position

%% Compute Bayesian recostruction over theta cycles in each running direction

% Split each cycle into 6 pieces
nPieces = 6; nBins = 200;
posAtCycle = interp1(positions(:,1),positions(:,2),cycles);
distanceTravelled = diff(posAtCycle,[],2);

% Reconstruct position separately for forward and backward direction:
cyclesCell = {cycles(distanceTravelled<0,:), cycles(distanceTravelled>0,:)}; % do forward cycles and backward cycles separately
titles = {'backward','forward'};
for i=1:2, 
    subplot(1,2,i);
    theseCycles = cyclesCell{i};
    splitCycles = SplitIntervals(theseCycles,'nPieces',nPieces);
    id = repmat((1:nPieces)',size(theseCycles,1),1); 
    % "id" this gives "ReconstructPosition" the location of the piece within a theta cycle. 
    % The function will group all the bins with "id==1" together in average(:,1). 
    % Plotting "average" therefore visualizes the average theta cycle.
    [estimations,actual,errors,average] = ReconstructPosition(positions,spikes,splitCycles,'training',behavior.run,'id',id,'nBins',nBins);
    % Here, "errors" is the "estimations" probability matrix visualized in 
    % the previous section, but centered on the actual position of the animal.
    yPosition = linspace(-1,1,size(estimations,1))';
    x = linspace(0,1,nPieces+1)'; x = mean([x(1:end-1) x(2:end)],2); x = [x; x+1];
    PlotColorMap(repmat(average,1,2),'y',yPosition,'x',x);
    set(gca,'xtick',0:1:2);
    xlabel('Theta cycles');
    PlotHVLines(0,'h','k--','linewidth',2);
    PlotHVLines(1,'v','k--','linewidth',2);
    clim([0 0.02])
    ylim([-1 1]*0.5);
    ylabel({'Decoded position relative to the true animal position','(underestimating -> overestimating)'});
    title(titles{i});
    set(gca,'box','off','TickDir','out','fontsize',12);
end

%% Compute Bayesian recostruction over all cycles

% You should see theta cycles sweeping towards future positions of the animal.
% This mimicks the running direction of the animal: late in the theta cycle,
% we see the position where the animal will be next, and when he's running 
% backward on the track (direction from 1 to 0), this is lower in the matrix.
% Note how if we hadn't restricted to "forward" or "backward" cycles and we had
% computed them together, a mix of these two kinds of cycles would have 
% resulted in a bad fit and theta cycles would not have been visible.
% To pool the two directions and plot them in the same matrix, we would
% need to flip the backward cycles so that positive positions are always in
% front of the animal.

nCycles = size(cycles,1); % use all running cycles
splitCycles = SplitIntervals(cycles,'nPieces',nPieces);
id = repmat((1:nPieces)',nCycles,1);
[estimations,actual,errors] = ReconstructPosition(positions,spikes,splitCycles,'training',behavior.run,'id',id,'nBins',nBins);
% "average" would not be useful here because around half the cycles need to be flipped first
errorArray = reshape(errors,size(errors,1),nPieces,[]); % make an array where the third dimension is cycleID
flip = distanceTravelled<0; % flip the matrices for backward running direction
errorArray(:,:,flip) = flipud(errorArray(:,:,flip)); % flip cycles with more reverse movement than forward
average = nanmean(errorArray,3);

clf
subplot(1,2,1);
PlotColorMap(average,'y',yPosition,'x',linspace(0,1,nPieces));
set(gca,'xtick',0:1);
xlabel('Theta cycles');
PlotHVLines(0,'h','k--','linewidth',1);
clim([0 0.02])
ylim([-1 1]*0.5);
ylabel({'Decoded position relative to the true animal position','(underestimating -> overestimating)'});
set(gca,'box','off','TickDir','out','fontsize',12);

% Find the best line using the Davidson method
% The method consists of trying out all possible linear fits to the probability
% matrix, keeping the fit that maximizes the sequence score. The score is the
% amount of probability within a specific distance threshold around the line.
% As Davidson, we will use a threshold of 15 bins away from the line, for a 
% 200 bin matrix, or 7.5%. Feel free to use adjust the threshold if your track
% is considerably different (threshold should be ~20 cm). 
threshold = ceil(0.075*nBins); 
sequence = ReplayScore(average,'threshold',threshold);
hold on; plot([0 1],yPosition([sequence.lineStart sequence.lineStop]),'w--','linewidth',2);
hold on; plot([0 1],yPosition([sequence.lineStart sequence.lineStop]-threshold),'w','linewidth',2);
hold on; plot([0 1],yPosition([sequence.lineStart sequence.lineStop]+threshold),'w','linewidth',2);

thisTitle = {['Average cycle: score=' num2str(sequence.score) ', p=' num2str(sequence.p) ', '],['weighted correlation = ' num2str(sequence.weightedCorrelation)]};
title(thisTitle);

subplot(1,2,2); cla
PlotColorMap(average,'y',yPosition,'x',x(1:nPieces));
set(gca,'xtick',0:1);
xlabel('Theta cycles');
PlotHVLines(0,'h','k--','linewidth',1);
clim([0 0.02])
ylim([-1 1]*0.5);
ylabel({'Decoded position relative to the true animal position','(underestimating -> overestimating)'});
set(gca,'box','off','TickDir','out','fontsize',12);
thisTitle = {'Average cycle:',['quadrant score =' num2str(sequence.quadrantScore)]};
plot([1/6 1/6],[-0.25 0],'k');
plot([0.5 0.5],[-0.25 0],'k');
plot([1/6 0.5],[0 0],'k');
plot([1/6 0.5],[-0.25 -0.25],'k');
plot([5/6 5/6],[0.25 0],'k');
plot([0.5 0.5],[0.25 0],'k');
plot([5/6 0.5],[0 0],'k');
plot([5/6 0.5],[0.25 0.25],'k');
plot([5/6 5/6],[-0.25 0],'w--');
plot([0.5 0.5],[-0.25 0],'w--');
plot([5/6 0.5],[0 0],'w--');
plot([5/6 0.5],[-0.25 -0.25],'w--');
plot([1/6 1/6],[0.25 0],'w--');
plot([0.5 0.5],[0.25 0],'w--');
plot([1/6 0.5],[0 0],'w--');
plot([1/6 0.5],[0.25 0.25],'w--');
title(thisTitle);
% The quadrant score is the mean probability within the black boxes
% (classic theta sequences, past-> future)
% minus the mean probability within the white boxes
% (which does not correspond to theta sequences and is a control)

%% Score individual theta cycles to get a distribution of scores

reconstructions = cell(nCycles,1);
for i=1:nCycles
    reconstructions{i} = errorArray(:,:,i);
end

sequences = ReplayScore(reconstructions{1},'threshold',threshold,'nShuffles',0); % no need for shuffles as we don't require a p-value for each individual theta cycle
for i=2:nCycles % append every theta sequence
    sequence = ReplayScore(reconstructions{i},'threshold',threshold,'nShuffles',0); % no need for shuffles (dramatically slower) as we don't require a p-value for each individual theta cycle
    sequences = MergeStructures(sequences,sequence);
end
% Fix the concatenation of nested structures:
jump.mean = cat(1,sequences.jump.mean);
jump.max = cat(1,sequences.jump.max);
sequences.jump = jump;
sequences = rmfield(sequences,'p'); % remove p-value as we did not compute it

% Plot statistics
clf
set(gcf,'position',[600 300 900 300]);
subplot(1,3,1);
histogram(sequences.quadrantScore,linspace(-1,1,50),'edgecolor','none');
PlotHVLines(0,'v','k--','linewidth',2);
p = signrank(sequences.quadrantScore); 
xlabel('quadrant score');
title(['Quadrant scores, p=' num2str(p)]);

subplot(1,3,2);
histogram(sequences.weightedCorrelation,linspace(-1,1,50),'edgecolor','none');
PlotHVLines(0,'v','k--','linewidth',2);
p = signrank(sequences.weightedCorrelation); 
title(['Weighted correlations, p=' num2str(p)]);
xlabel('weighted correlation');

subplot(1,3,3);
histogram(sequences.slope,linspace(-1,1,50)*nBins/nPieces,'edgecolor','none');
PlotHVLines(0,'v','k--','linewidth',2);
p = signrank(sequences.slope); 
title(['Slopes, p=' num2str(p)]);
xlabel('slope (spatial bins / cycle)');

%% Save theta sequences

% Before saving the theta sequences, I recommend you save some additional
% details which can be useful:

sequences.reconstructions = reconstructions;
sequences.average = average;
sequences.intervals = cycles;
% Add "state" as a standard measure that will also appear in replay
% sequence structures. "1" is code for the awake state (as seen in "stateNames")
% Since cycles are restricted to running, all theta sequences are awake data
sequences.state = ones(nCycles,1); 
sequences.stateNames = {'WAKE','NREM','REM'};
% The subsession within "MergePoint" that the cycles belong to:
try
    MergePoints = getStruct(basepath,'MergePoints');
    [~,epoch] = InIntervals(cycles,MergePoints.timestamps);
    sequences.epoch = epoch;
    sequences.epochNames = MergePoints.foldernames';
end
% More details about the events:
sequences.events.timestamps = cycles;
sequences.events.peaks = peaks;
sequences.events.amplitude = amplitude;
try
    rt = RelativeTimes(spikes(:,1),cycles);
    [m,r,p] = CircularMean(rt*2*pi,1,spikes(:,2));
    sequences.phaseLocking.angle = m; sequences.phaseLocking.meanVectorLength = r; sequences.phaseLocking.p = p;
end
sequences.params.binsPerCycle = nPieces;
sequences.params.project = 'OML'; % feel free to label your own project
sequences.params.environment = 'LinearTrack';
sequences.training.intervals = behavior.run; % the training intervals that were used to decode the position
sequences.training.positions = positions; % the positions that were used to decode the position

% Save sequences using the following names:
thetaSequences = sequences;
thetaFilename = fullfile(basepath,[basename '.thetaSequences.mat']);
save(thetaFilename,'thetaSequences');


