% This function is to help with the tutorial pipelineBatch

function [mPETH1,mPETH2,basepath] = pipelineBatch_BatchExample(basepath,varargin)

% batch = StartBatch(@BatchExample,'OMLlinear.batch'); % This is a little comment that helps me remember how to call this function

% This is a simple function which would load some basic data and output the mean response of each unit to 
% pre-task sleep ripples and post-task sleep ripples

% Note that "basepath" is also an output of this function. I find it useful to make sure I never get lost and
% I can always find out which session particular results come from.


disp([datestr(clock) ': Starting session ' basepath '...']);
basename = basenameFromBasepath(basepath);
ripples = getStruct(basepath,'ripples');
MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
% Load NREM sleep epochs:
load(fullfile(basepath,[basename '.SleepState.states.mat']),'SleepState');
NREM = bsxfun(@plus,SleepState.ints.NREMstate,[-1 0]);
spikeStructure = importSpikes('basepath',basepath,'CellType',"Pyramidal Cell",'brainRegion','CA1'); % load CA1 pyramidal cells
nUnits = length(spikeStructure.times);
preSleep = OverlapIntervals(NREM,sleep(1,:));
preSleep = TruncateIntervals(preSleep,3600);
postSleep = OverlapIntervals(NREM,sleep(2:end,:));
postSleep = TruncateIntervals(postSleep,3600);

preRipples = Restrict(ripples.timestamps(:,1),preSleep);
postRipples = Restrict(ripples.timestamps(:,1),postSleep);

% For each unit, compute the mean PETH in pre-task sleep and post-task sleep
nBins = 101; durations = [-1 1];
[mPETH1,mPETH2] = deal(nan(nUnits,nBins));
for i=1:nUnits
    s = spikeStructure.times{i};
    mPETH1(i,:) = mPETH(s,preRipples,'durations',durations,'nBins',nBins);
    mPETH2(i,:) = mPETH(s,postRipples,'durations',durations,'nBins',nBins);
end









