% script_OML_fixCheeseBehaviorFile
function basepath = script_OML_fixCheeseBehaviorFile(basepath,condition,taskType)
%% Basic pipeline for place field analysis: adapted for the OML project (custom linearization, 100 bins, smooth 1, threshold = 0.2)
% These are the basic fucntions for performing place field analysis (in
% 1D for now). Is still work in progress and anyone is welcome to contribute.

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
% 
% try
%     load(fullfile(basepath,[basename '.placeFields.cellinfo.mat']),'placeFieldStats');
%     disp(['Place fields for ' basepath ' are already computed. Aborting...']);
% %     keyboard
%     return
% end



%%


if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath); 

load(fullfile(basepath,[basename '_hmmLinearized.mat']));
behavior = getStruct(basepath,'animal.behavior');
l = hmmLinearized.linearized;
behavior.position.linearized = l(:)';

disp([datestr(clock) ': Starting session ' basepath]);

cd(basepath);

% save(fullfile(basepath,[basename '.animalBACKUP.behaviorBACKUP.mat']),'behavior');
MergePoints = getStruct(basepath,'MergePoints');

cell_metrics = getStruct(basepath,'cell_metrics');
% spikesStruct = getStruct(basepath,'spikes.cellinfo');
spikesCell = cell_metrics.spikes.times';
maze = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'linear')),MergePoints.foldernames),:));
if isempty(maze)
    maze = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'run')),MergePoints.foldernames),:));
end
nSpikes = cellfun(@(x) sum(InIntervals(x,maze)), spikesCell);
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
pyr = pyr & nSpikes>200;
pyrID = find(pyr);

spikesCell = cell_metrics.spikes.times;

try
    error
    run = behavior.run;
catch
    t = (behavior.timestamps(1):mode(diff(behavior.timestamps(:))):behavior.timestamps(end))';
    % Old way to define running:
    if strcmp(projectName,'OML18') % there is an issue with behavior.position.x in OML18
        interpolated = sortrows([behavior.positionTrials{1};behavior.positionTrials{2}(:,1) 1-behavior.positionTrials{2}(:,2)]);
        interpolated = nani(interpolated); interpolated(:,3) = interpolated(:,2);
        interpolated(diff(interpolated(:,1))==0,:) = [];
    else
        ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:));
        interpolated = interp1(behavior.timestamps(ok)',[behavior.timestamps(ok)' behavior.position.x(ok)',behavior.position.y(ok)'],t);
        interpolated(isnan(interpolated(:,1)),:) = [];
    end
    speed = LinearVelocity(interpolated,5); t = speed(:,1);
    try
    allPositions = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}]);
    catch
        allPositions = [behavior.timestamps(:) behavior.position.linearized(:)]; allPositions(:,2) = ZeroToOne(allPositions(:,2));
    end
    [up,down] = ZeroCrossings([allPositions(:,1),allPositions(:,2)-0.5]);
    midpoints = allPositions(up|down,1); 
    ok = InIntervals(allPositions(FindClosest(allPositions(:,1),midpoints),2),[0.3 0.7]); % true midpoints are close to positions around 0.5 (not 1 to 0 wrap crossings)
    midpoints = midpoints(ok);
    topSpeed = interp1(speed(:,1),speed(:,2),midpoints);
    threshold = nanmedian(topSpeed)/10; % the threshold is 10% of the peak speed
    run = t(FindInterval(speed(:,2)>threshold));
    run = ConsolidateIntervals(run,'epsilon',0.01);

    % Make sure periods without position data are not erroneously included in "running"
    t0 = (behavior.timestamps(1):mode(diff(behavior.timestamps(:))):behavior.timestamps(end))';
    gaps = abs(reshape(behavior.timestamps(FindClosest(behavior.timestamps(:),t0(:))),[],1) - t0(:));
    noData = t0(FindInterval(gaps>0.1)); if size(noData,2)==1, noData=noData'; end
    run = SubtractIntervals(run,noData);

    % Make sure our each run epoch is confined within the same trial: no run epoch should begin in one trial and continue over the next trial
    run = SubtractIntervals(run,bsxfun(@plus,sort(behavior.trials(:)),[-1 1]*0.00001));

    [in,w] = InIntervals(t,run);
    peak = Accumulate(w(in),speed(in,2)','mode','max');

    % remove run epochs that don't reach the speed threshold
    if strcmp(projectName,'OLM21') || strcmp(projectName,'OML22') || strcmp(projectName,'OML23')  % there is an issue with behavior.position.x in OML18
        run(peak<20,:) = [];
    elseif strcmp(projectName,'OML18') || strcmp(projectName,'OML19')
        run(peak<0.1,:) = [];
    else
        run(peak<20,:) = [];
    end

    MergePoints = getStruct(basepath,'MergePoints');
    sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
    run(IntervalsIntersect(run,sleep),:) = [];
    run(diff(run,[],2)<0.5,:) = []; % remove run epochs lasting for less than 0.5s

    behavior.run = run;
end

% posTrials = behavior.positionTrials;

load('cheesbTrialsDayOffset.mat');
behavior.trials = [cheesbTrialsDay.trials.start cheesbTrialsDay.trials.end];
trialKind = 2 - ones(size(behavior.trials,1),1)*condition;
behavior.trialID(:,1) = trialKind;
if condition==1, stim = [behavior.trials(1,1) behavior.trials(end,end)]; else, stim = zeros(0,2); end
behavior.stimON = stim;

%%
j=1;
posTrials{j} = Restrict([behavior.timestamps(:),behavior.position.linearized(:)],behavior.trials); 

[in,w] = InIntervals(posTrials{j},behavior.trials);
ok = find(in & InIntervals(posTrials{j},run));
[mins,that] = Accumulate(w(ok),posTrials{j}(ok,2),'mode','min'); % min position for each lap
maxes = Accumulate(w(ok),posTrials{j}(ok,2),'mode','max'); % max position for each lap
mm = [nanmedian(mins) nanmedian(maxes)]; % define starting position (0) and end position (1)
new = (posTrials{j}(:,2)-mm(1))./diff(mm);
new(new<0) = 0; new(new>1) = 1;
posTrials{j}(:,2) = new;
posTrials{j} = posTrials{j}(in & ~isnan(posTrials{j}(:,2)),:);
this = Restrict(posTrials{j},run);
% make sure positions are going from 0 to 1, not from 1 to 0!
d = diff(this(:,2));
if nanmedian(d(abs(d)>0))<0
    posTrials{1}(:,2) = 1-posTrials{j}(:,2);
end

behavior.positionTrials = posTrials;

%%

posTrials{j} = Restrict(posTrials{j},run);
% cut away the points with saturated positions:
saturated = find(posTrials{j}(:,2)==0);
bad = saturated(diff(saturated)==1); % consecutive saturated points will be removed (keep last point only)
posTrials{j}(bad,:) = [];
saturated = find(posTrials{j}(:,2)==1);
bad = saturated([false;diff(saturated)==1]); % consecutive saturated points will be removed (keep first point only)
posTrials{j}(bad,:) = [];

behavior.positionTrialsRun = posTrials;
save(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');
