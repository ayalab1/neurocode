function [g,hs,basepath] = BatchV1ProbeTrials(basepath)

%%
rng(1);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);

session = getStruct(basepath,'session');
ripples = getStruct(basepath,'ripples');
MergePoints = getStruct(basepath,'MergePoints');
EMG = getStruct(basepath,'EMG');
immobility = EMG.timestamps(FindInterval(EMG.data<0.6)); immobility(diff(immobility,[],2)<1,:) = [];
postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));

[spikes,regionID,regionNames,spikesCell] = GetAyaSpikes(basepath);
% Get a single string to include in figure labels
regionListString = regionNames; regionListString(2,:) = repmat({'->'},1,size(regionListString,2)); regionListString{end} = '';
regionListString = strcat(regionListString{:});

regionNames = regionNames(:); regionID = regionID(:);
regionsCell = regionNames(regionID);

%%
clf
matrix = dlmread('stimuli_synced.times');
matrix(isnan(matrix(:,1)),:) = [];

nUnits = length(spikesCell);
% stim = matrix(diff([0;matrix(:,2)])~=0 & matrix(:,2)>-1,:);
ok = diff([0;matrix(:,2)])~=0 & matrix(:,end)>0 & matrix(:,3)==0 & matrix(:,2)>-2;
stims = matrix(ok,:);
u = unique(stims(:,2)); %u(u==-2) = [];


%%
ok = diff([0;matrix(:,2)])~=0 & matrix(:,end)>0 & matrix(:,3)==0 & matrix(:,2)>-2; stims = matrix(ok,:);
upper = Accumulate(stims(:,end),stims(:,2),'mode','max');
audiovisual = find(upper==90);
if length(audiovisual)>1, n = Accumulate(stims(:,end)); audiovisual = audiovisual(findmax(n(audiovisual))); end
q = stims(ismember(stims(:,end),audiovisual),:);
if any(q(:,2)<0), error('blank trials present'); end
q(q(:,2)<0,:) = [];
nRepetitions = 30; nStimuli = 2; stimuli = [0 90]; nProbeTrials = 5; firstPossibleProbeTrial = 10; 
q(nRepetitions*nStimuli+1:end,:) = [];
trialType = q(:,2); 
for kk = 1:1000,
    % Select random trials to be probe trials
    probeTrial = false(size(trialType));
    if nProbeTrials>0
        %     rng(floor(datenum(clock))); % set a new random number generator every day. NOTE: this was added on May 05 2022, and it was therefore not applicable to Jean and AO52
        allPossibleLines = 1:nRepetitions*nStimuli;
        % write "-1" in 3rd column for probe trials:
        ok = false;
        while ~ok % start over if you got stuck in an infinite loop
            probeTrial = false(size(trialType));
            try
            for i = 1:length(stimuli)
                if ~(stimuli(i)>-1) continue; end % there is no sound before the control stimulus anyway (or the blank screen (NaN))
                possibleLines = allPossibleLines(trialType==stimuli(i));
                possibleLines(1:firstPossibleProbeTrial-1) = []; % first X trials are off limits and should not be probe trials
                ok = false;
                timestamp = GetSecs;
                while ~ok && GetSecs - timestamp<2 % if you haven't found a solution in 2s, then we're probably stuck and should start over
                    indices = sort(randperm(length(possibleLines),nProbeTrials));
                    lines = possibleLines(indices)';
                    ok = true;
                    if any(diff(indices) == 1), ok = false; end % avoid having 2 probe trials in a row (for the same stimulus)
                    if any(probeTrial(ismember(allPossibleLines,[lines-1; lines+1]))), ok = false; end % avoid having 2 probe trials in a row (even for different stimuli)
                end
                probeTrial(lines) = -1; % selected for a probe trial
            end
            end
        end
    end
    probeTrials(kk,:) = probeTrial';
end

[h0,ht] = PETH(spikes(:,1),q(:,1)-5,'durations',[-1 1]*0.5,'nBins',201);
h0 = h0/mode(diff(ht));
response0 = mean(h0(:,InIntervals(ht,[0.09 0.14])),2);
[h,ht] = PETH(spikes(:,1),q(:,1),'durations',[-1 1]*0.5,'nBins',201);
h = h/mode(diff(ht));
response = mean(h(:,InIntervals(ht,[0.09 0.14])),2);

this = probeTrials*response;
worstCase = probeTrials(findmin(this),:)';

g = Group(response(1:firstPossibleProbeTrial*2), response(worstCase),response0);
hs = Group(h(1:firstPossibleProbeTrial*2,:), h(worstCase,:),h0);




