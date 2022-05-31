% script_OML_pipelinePlaceCells
function script_OML_pipelinePlaceCells(basepath)
%% Basic pipeline for place field analysis: adapted for the OML project (custom linearization, 100 bins, smooth 1, threshold = 0.2)
% These are the basic fucntions for performing place field analysis (in
% 1D for now). Is still work in progress and anyone is welcome to contribute.

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

try
    load(fullfile(basepath,[basename '.placeFields.cellinfo.mat']),'placeFieldStats');
    disp(['Place fields for ' basepath ' are already computed. Aborting...']);
%     keyboard
%     return
end

%%

disp([datestr(clock) ': Starting session ' basepath]);

cd(basepath);

behavior = getStruct(basepath,'animal.behavior');

cell_metrics = getStruct(basepath,'cell_metrics');
% spikesStruct = getStruct(basepath,'spikes.cellinfo');
spikesCell = cell_metrics.spikes.times';
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
    else
        ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:));
        interpolated = interp1(behavior.timestamps(ok)',[behavior.timestamps(ok)' behavior.position.x(ok)',behavior.position.y(ok)'],t);
        interpolated(isnan(interpolated(:,1)),:) = [];
    end
    speed = LinearVelocity(interpolated,5); t = speed(:,1);
    allPositions = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}]);
    [up,down] = ZeroCrossings([allPositions(:,1),allPositions(:,2)-0.5]);
    midpoints = allPositions(up|down,1);
    topSpeed = interp1(speed(:,1),speed(:,2),midpoints);
    threshold = nanmedian(topSpeed)/10; % the threshold is 10% of the peak speed
    run = t(FindInterval(speed(:,2)>threshold));
    run = ConsolidateIntervals(run,'epsilon',0.01);

    % Make sure periods without position data are not erroneously included in "running"
    t0 = (behavior.timestamps(1):mode(diff(behavior.timestamps(:))):behavior.timestamps(end))';
    gaps = abs(reshape(behavior.timestamps(FindClosest(behavior.timestamps(:),t0(:))),[],1) - t0(:));
    noData = t0(FindInterval(gaps>0.1));
    run = SubtractIntervals(run,noData);

    % Make sure our each run epoch is confined within the same trial: no run epoch should begin in one trial and continue over the next trial
    run = SubtractIntervals(run,bsxfun(@plus,sort(behavior.trials(:)),[-1 1]*0.00001));

    [in,w] = InIntervals(t,run);
    peak = Accumulate(w(in),speed(in,2)','mode','max');

    % remove run epochs that don't reach the speed threshold
    if strcmp(projectName,'OLM21') || strcmp(projectName,'OML22') % there is an issue with behavior.position.x in OML18
        run(peak<20,:) = [];
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
if size(behavior.trialID,2)>1, trialKind = 2-behavior.trialID(:,2); % trialKind=1 is for stim trials, trialKind=1 is for Off trials
    1
else trialKind = behavior.trialID(:,1); 2
    if max(trialKind)==1, trialKind = 2-trialKind; 3
    end
end 
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
    mm = [median(mins) median(maxes)]; % define starting position (0) and end position (1)
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
intervalsCell = {stimIntervals,nonstimIntervals};
nUnits = length(spikesCell);
for j=1:2 % Restrict to running epochs
    intervals = intervalsCell{j};
    restrictedPos = Restrict(posTrials{j},intervals,'shift','on'); 
    restrictedSpiketimes = cell(1,nUnits); for i=1:nUnits, restrictedSpiketimes{i} = Restrict(spikesCell{i},intervals,'shift','on'); end
    restrictedSpikestruct.times = restrictedSpiketimes; restrictedSpikestruct.UID = 1:nUnits;
    if j==1, firingMaps = firingMapAvg({restrictedPos},restrictedSpikestruct,'nBins',100,'smooth',1,'savemat',false);
    else
        temp = firingMapAvg({restrictedPos},restrictedSpikestruct,'nBins',100,'smooth',1,'savemat',false);
        for i=1:nUnits, firingMaps.rateMaps{i}{j} = temp.rateMaps{i}{1}; firingMaps.countMaps{i}{j} = temp.countMaps{i}{1}; firingMaps.occupancy{i}{j} = temp.rateMaps{i}{1}; end
    end
end
save(fullfile(basepath,[basename '.firingMapsAvg.cellinfo.mat']),'firingMaps');

%% 1.1- Save struct for cell explorer (you can display ratemaps in cell explorer)
% Optimized for only 1d, though 2d works too
ratemap = firingMaps;
for i = 1:length(ratemap.rateMaps)
    ratemap.map{i} = ratemap.rateMaps{i}{1};
end
ratemap.x_bins = ratemap.params.x;
ratemap.y_bins = ratemap.params.y;
ratemap.rateMaps = [];
save(fullfile(basepath,[basename '.ratemap.firingRateMap.mat']),'ratemap');

%% 2- Detect (average) place fields
[placeFieldStats] = findPlaceFieldsAvg1D('firingMaps',firingMaps,'minPeak',1,'sepEdge',0.04,'threshold',0.2,'doPlot',false);
% save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');

% % 2.1 alternative to detect place fields 2D
% for i = 1:length(firingMaps.rateMaps)
%     map.z = firingMaps.rateMaps{i}{1};
%     map.time = firingMaps.occupancy{i}{1};
%     map.count = firingMaps.countMaps{i}{1};
%     map.x = firingMaps.params.x;
%     map.y = firingMaps.params.y;
%     placeFieldStats{i} = MapStats(map);
% end
% save(fullfile(basepath,[basename '.placeFields.cellinfo.mat']),'placeFieldStats');

% %% 3- Calculate template of place fields in the maze (for decoding, etc.)
% [placeFieldTemplate] = findPlaceFieldsTemplate('placeFieldStats',placeFieldStats,'firingMaps',firingMaps);
% % save([basename '.placeFieldTemplate.mat'],'placeFieldTemplate');

% %% 4- Phase precession
% % theta phase
% 
% [parentFolder,dayName] = fileparts(basepath);
% [~,animalName] = fileparts(parentFolder);
% switch animalName
%     case 'OML18'
%         channel = 34;
%     case 'OML19'
%         channel = 34;
%     case 'OLM21'
%         channel = 22;
%     case 'OML22'
%         channel = 163;
%     otherwise
%         keyboard
% end
% 
% lfp = getLFP(channel+1);
% theta = bz_Filter(lfp,'passband',[5 15]);
% 
% % boundaries of each PF
% for i=1:numel(spikesStruct.UID) %
%     for j=1:2
%         for k=1:length(placeFieldStats.mapStats{i}{j}.peak)
%             if placeFieldStats.mapStats{i}{j}.peak(k) ~= 0
%                 boundaries{i}{j}(k,1)= firingMaps.params.x(placeFieldStats.mapStats{i}{j}.fieldX(k,1));
%                 boundaries{i}{j}(k,2)= firingMaps.params.x(placeFieldStats.mapStats{i}{j}.fieldX(k,2));
%             else
%                 boundaries{i}{j}(k,1)= NaN;
%                 boundaries{i}{j}(k,2)= NaN;
%             end
%         end
%     end
% end
% 
% count = 0;
% phases = [theta.timestamps, theta.phase];
% % calculate phase precession
% tic
% intervalsCell = {stimIntervals,nonstimIntervals}; intervals = intervalsCell{j};
% for i=1:numel(spikesStruct.UID)
%     for j=1:2
%         for k=1:length(placeFieldStats.mapStats{i}{j}.x) % number of place fields
%             if ~isnan(placeFieldStats.mapStats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed
%                 count=count+1;
%                 this = posTrials{j}; mm = [min(this(:,2)) max(this(:,2))];
%                 this(:,2) = (this(:,2)-mm(1))./diff(mm);
%                 theseBoundaries = (boundaries{i}{j}(k,:)-mm(1))./diff(mm); theseBoundaries(theseBoundaries<0) = 0; theseBoundaries(theseBoundaries>1) = 1;
%                 if any(isnan(theseBoundaries)), continue; end
%                 if numel(Restrict(spikesStruct.times{i},intervals))<20, continue; end
%                 intervals = intervalsCell{j};
%                 [dataPP{i,j},statsPP{i,j}] = PhasePrecession(this,Restrict(spikesStruct.times{i},intervals),phases,'boundaries',theseBoundaries); 
%             end
%         end
%     end
%     disp(['Done with unit ' num2str(i) '/' num2str(numel(spikesStruct.UID)) '; estimated time remaining=' num2str(toc/i*(numel(spikesStruct.UID)-i))]);
% end
% 
% PP = struct; PP.dataPP = dataPP; PP.statsPP = statsPP;
% save(fullfile(basepath,[basename '.phasePrecession.mat']),'PP');
% 
% % for i=1:numel(dataPP)
% %     if ~isempty(dataPP{i})
% %      PlotPhasePrecession(dataPP{i},statsPP{i});
% %      pause;
% %      close all;
% %     end
% % end
% 
% nBins = [10 25];
% for i=1:size(dataPP,1)
%     for j=1:size(dataPP,2)
%         if isempty(dataPP{i,j}), continue; end
%         qq = [(dataPP{i,j}.position.x - statsPP{i,j}.boundaries(1))./diff(statsPP{i,j}.boundaries) dataPP{i,j}.position.phase];
%         this = Restrict(qq,[0 1]); this(:,2) = wrap(this(:,2),2);
%         this(:,1) = ceil(this(:,1)*nBins(1));
%         this(:,2) = ceil(this(:,2)*nBins(2)/(2*pi)); this(this==0) = 1;
%         saved{i,j} = Accumulate(this,1,'size',nBins)';
%         smoothed{i,j} = Smooth(saved{i,j}/sum(saved{i,j}(:)),1,'type','cc');
%     end
% end
% 
% figure(10); clf
% x = linspace(0,1,nBins(1));
% y = linspace(0,4*pi,nBins(2)*2)/pi;
% for j=1:2
%     q = cat(3,smoothed{:,j});
%     subplot(1,2,j); PlotColorMap(repmat(nanmean(q,3),2,1),'x',x,'y',y);
% end
% drawnow

