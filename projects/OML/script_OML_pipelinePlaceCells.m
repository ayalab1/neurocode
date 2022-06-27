% script_OML_pipelinePlaceCells
function basepath = script_OML_pipelinePlaceCells(basepath)
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

disp([datestr(clock) ': Starting session ' basepath]);

cd(basepath);

behavior = getStruct(basepath,'animal.behavior');
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
    stim = behavior.stimON;
catch
    % convert from old formats (for OML18 and OML19)
    load(fullfile(basepath,'stimInt.mat'), 'stimInt');
    stim = stimInt.On;
    behavior.stimON = stim;
    try
        load(fullfile(basepath,[basename '_Pos.mat']), 'Pos');
        takeColumn = (diff(sum(Pos.ON(:,2:3)>0),[],2)>0)+2;
        posOn = Pos.ON(:,[1 takeColumn]);
        posOn(:,2) = posOn(:,2)-min(posOn(:,2)); posOn(:,2) = posOn(:,2)/max(posOn(:,2)); % normalize (separately as this was badly done for some directions)
        takeColumn = (diff(sum(Pos.Off(:,2:3)>0),[],2)>0)+2;
        posOff = Pos.Off(:,[1 takeColumn]);
        posOff(:,2) = posOff(:,2)-min(posOff(:,2)); posOff(:,2) = posOff(:,2)/max(posOff(:,2)); % normalize (separately as this was badly done for some directions)
        posOff(:,2) = 1-posOff(:,2); % flip because the code below expects it
    catch  % OML18 days 5 6 7 don't have either structure, so load posTrials in those...
        load(fullfile(basepath,['posTrials.mat']), 'posTrials');
        posOn = posTrials{2}; posOff = posTrials{1};
    end
    q = sortrows([posOn; posOff]);
    behavior.position.x = q(:,2)';
    behavior.position.y = ones(size(q(:,2)))';
    behavior.positionTrials{1} = posOn;
    behavior.positionTrials{2} = posOff;
    isBoundary = zscore(diff(posOn(:,1)))>1;
    trialsOn = [posOn([true; isBoundary],1) posOn([isBoundary;true],1)];
    isBoundary = zscore(diff(posOff(:,1)))>1;
    trialsOff = [posOff([true; isBoundary],1) posOff([isBoundary;true],1)];
    trials = Group(trialsOn,trialsOff);
    trials(:,3) = 2-trials(:,3);
    behavior.trials = trials(:,1:2);
    behavior.trialID = trials(:,[3 3]);
end

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
        keyboard
    end

    MergePoints = getStruct(basepath,'MergePoints');
    sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
    run(IntervalsIntersect(run,sleep),:) = [];
    run(diff(run,[],2)<0.5,:) = []; % remove run epochs lasting for less than 0.5s

    behavior.run = run;
end

% posTrials = behavior.positionTrials;

% normalize positions from 0 to 1, with 0 and 1 corresponding to the actual beginning/end of a lap and not an outlier point
if size(behavior.trialID,2)>1, trialKind = 2-behavior.trialID(:,2); % trialKind=1 is for stim trials, trialKind=0 is for Off trials
else trialKind = behavior.trialID(:,1);
    if max(trialKind)==1, trialKind = 2-trialKind;
    end
end 
behavior.trialID(:,1) = trialKind;
nonrun = SubtractIntervals([0 Inf],run); nonstim = SubtractIntervals([0 Inf],stim);
stimIntervals = SubtractIntervals(behavior.trials(trialKind==1,:),ConsolidateIntervals(sortrows([nonrun; nonstim]))); % only stim run periods of the correct trial type
nonstimIntervals = SubtractIntervals(behavior.trials(trialKind==2,:),ConsolidateIntervals(sortrows([nonrun; stim])));  % only stim run periods

%%
for j=1:2,posTrials{j} = Restrict([behavior.timestamps(:),behavior.position.linearized(:)],behavior.trials(trialKind==j,:)); end

% check if we need to flip "posTrials":
if sum(InIntervals(Restrict(posTrials{1},run),stim)) < sum(InIntervals(Restrict(posTrials{2},run),stim)), posTrials = posTrials([2 1]); disp('flipped'); end
for j=1:length(posTrials),
    [in,w] = InIntervals(posTrials{j},behavior.trials(trialKind==j,:));
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
        posTrials{j}(:,2) = 1-posTrials{j}(:,2);
    end
end
behavior.positionTrials = posTrials;

%%

for j=1:2,
    posTrials{j} = Restrict(posTrials{j},run); 
    % cut away the points with saturated positions:
    saturated = find(posTrials{j}(:,2)==0);
    bad = saturated(diff(saturated)==1); % consecutive saturated points will be removed (keep last point only)
    posTrials{j}(bad,:) = [];
    saturated = find(posTrials{j}(:,2)==1);
    bad = saturated([false;diff(saturated)==1]); % consecutive saturated points will be removed (keep first point only)
    posTrials{j}(bad,:) = [];
end

behavior.positionTrialsRun = posTrials;
save(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');

%% 1-Calculate firing maps

posTrials = behavior.positionTrialsRun;

% intervalsCell = {stimIntervals,nonstimIntervals};

nonrun = SubtractIntervals([0 Inf],run); nonstim = SubtractIntervals([0 Inf],stim);
% for j=1:2, intervalsCell{j} = SubtractIntervals(behavior.trials(trialKind==j,:),ConsolidateIntervals(sortrows([nonrun]))); end
nUnits = length(spikesCell);
nBins = 100;
clear firingMaps placeFieldStats
firingMaps.UID = 1:nUnits;
placeFieldStats.UID = 1:nUnits;
firingMaps.sessionName = sessionID;
firingMaps.params.mode = 'discard'; firingMaps.params.nBins = nBins; 
for j=1:2 % Restrict to running epochs
    %     intervals = intervalsCell{j};
    %     restrictedPos = Restrict(posTrials{j},intervals,'shift','on');
    %     restrictedSpiketimes = cell(1,nUnits);
    for i=1:nUnits
        ii = i;
        %         restrictedSpiketimes{i} = Restrict(spikesCell{i},intervals,'shift','on');

        intervals = SubtractIntervals(behavior.trials(trialKind==j,:),ConsolidateIntervals(sortrows([nonrun])));
        restrictedPos = Restrict(posTrials{j},intervals,'shift','on');
        s = Restrict(spikesCell{ii},intervals,'shift','on');
        [curve,stats] = FiringCurve(restrictedPos,s,'nBins',nBins);
        firingMaps.rateMaps{i,1}{j} = curve.rate;
        firingMaps.countMaps{i,1}{j} = curve.count;
        firingMaps.occupancy{i,1}{j} = curve.time;
        confirmed = false(size(stats.x));
        % check if place field is preserved in each half of the data
        if ~isnan(stats.x(1)) || pyr(i)==false % do this only if a place field was found
            fields = double(permute(stats.field,[2 3 1]));
            firingMaps.params.x = curve.x; firingMaps.params.y = curve.y;

            half = find(trialKind==j); half = half(1:2:end); intervals = SubtractIntervals(behavior.trials(half,:),ConsolidateIntervals(sortrows([nonrun])));
            restrictedPos = Restrict(posTrials{j},intervals,'shift','on');
            s = Restrict(spikesCell{ii},intervals,'shift','on');
            [~,stats1] = FiringCurve(restrictedPos,s,'nBins',nBins);

            half = find(trialKind==j); half = half(2:2:end); intervals = SubtractIntervals(behavior.trials(half,:),ConsolidateIntervals(sortrows([nonrun])));
            restrictedPos = Restrict(posTrials{j},intervals,'shift','on');
            s = Restrict(spikesCell{ii},intervals,'shift','on');
            [~,stats2] = FiringCurve(restrictedPos,s,'nBins',nBins);

            if ~isnan(stats1.x(1)) && ~isnan(stats2.x(1)), % if no field was found, none of the fields are confirmed
                fields1 = double(permute(stats1.field,[2 3 1]));
                [c1,p1] = nancorr(fields,fields1);
                confirmed = any(c1>0 & p1<0.05,2);
                fields2 = double(permute(stats2.field,[2 3 1]));
                [c2,p2] = nancorr(fields,fields2);
                confirmed = confirmed & any(c2>0 & p2<0.05,2);
            end
        end
        stats.field = permute(stats.field,[2 3 1]);
        if ~any(confirmed),
            stats.field = false(nBins,0); stats.x = nan; stats.y = nan; stats.fieldX = nan(1,2); stats.fieldY = nan(1,2);
            stats.size = 0; stats.peak = 0; stats.mean = 0;
        else
            stats.field = stats.field(:,confirmed);  stats.x = stats.x(confirmed); stats.y = stats.y(confirmed); stats.fieldX = stats.fieldX(confirmed,:); stats.fieldY = stats.fieldY(confirmed,:);
            stats.size = stats.size(confirmed);  stats.peak = stats.peak(confirmed);  stats.mean = stats.mean(confirmed);
        end
        placeFieldStats.mapStats{i,1}{j} = stats;

        curves{j}(i,:) = curve.rate;
        %         restrictedSpikestruct.times = restrictedSpiketimes; restrictedSpikestruct.UID = 1:nUnits;
        %         if j==1, firingMaps = firingMapAvg({restrictedPos},restrictedSpikestruct,'nBins',100,'smooth',1,'savemat',false);
        %         else
        %             temp = firingMapAvg({restrictedPos},restrictedSpikestruct,'nBins',100,'smooth',1,'savemat',false);
        %             for i=1:nUnits, firingMaps.rateMaps{i}{j} = temp.rateMaps{i}{1}; firingMaps.countMaps{i}{j} = temp.countMaps{i}{1}; firingMaps.occupancy{i}{j} = temp.rateMaps{i}{1}; end
    end
end
save(fullfile(basepath,[basename '.firingMapsAvg.cellinfo.mat']),'firingMaps');
save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');

%% Export for cell explorer (you can display ratemaps in cell explorer)
% Optimized for only 1d, though 2d works too
ratemap = firingMaps;
for i = 1:length(ratemap.rateMaps)
    ratemap.map{i} = ratemap.rateMaps{i}{1};
end
ratemap.x_bins = ratemap.params.x;
ratemap.y_bins = ratemap.params.y;
ratemap.rateMaps = [];
save(fullfile(basepath,[basename '.ratemap.firingRateMap.mat']),'ratemap');
