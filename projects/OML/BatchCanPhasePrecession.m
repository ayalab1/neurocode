function [specStim,specNonStim,f,pyr,basepath] = BatchCanPhasePrecession(basepath)

% /home/programs/Matlab/train/SessionsReactivations1-Channels.batch
% RESULT: as expected, RT rather than filtered phase time changes nothing
%%

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,animalName] = fileparts(parentFolder);
sessionID = [animalName '_' dayName];

tic
name = [sessionID '_ThetaSpectrum.mat'];
try
    load(['M:\home\raly\Documents\code\new\tempFiles\' name],'specStim','specNonStim','f');
    f = linspace(0,3,400)';
    raw = specStim; 
    specStim = interp1(linspace(0,3,size(raw,2))',raw',f)'; 
    raw = specNonStim; 
    specNonStim = interp1(linspace(0,3,size(raw,2))',raw',f)'; 

    cell_metrics = getStruct(basepath,'cell_metrics');
    pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';

    return
end

%%

cd(basepath)
basepath
% thetaPhase = getStruct(basepath,'theta','phases');
% ok = ~isnan(thetaPhase(:,2)); thetaPhase(ok,2) = unwrap(thetaPhase(ok,2));

[spikes,regionID,regionNames,spikesCell,order] = GetAyaSpikes(basepath);
behavior = getStruct(basepath,'animal.behavior');
MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
try
    stim = behavior.stimON;
catch
    % convert from old formats (for OML18 and OML19)
    load(fullfile(basepath,'stimInt.mat'), 'stimInt'); stim = stimInt.On;
    try
        load(fullfile(basepath,[basename '_Pos.mat']), 'Pos');
        takeColumn = (diff(sum(Pos.ON(:,2:3)>0),[],2)>0)+2; posOn = Pos.ON(:,[1 takeColumn]);
        posOn(:,2) = posOn(:,2)-min(posOn(:,2)); posOn(:,2) = posOn(:,2)/max(posOn(:,2)); % normalize (separately as this was badly done for some directions)
        takeColumn = (diff(sum(Pos.Off(:,2:3)>0),[],2)>0)+2; posOff = Pos.Off(:,[1 takeColumn]);
        posOff(:,2) = posOff(:,2)-min(posOff(:,2)); posOff(:,2) = posOff(:,2)/max(posOff(:,2)); % normalize (separately as this was badly done for some directions)
        posOff(:,2) = 1-posOff(:,2); % flip because the code below expects it
    catch  % OML18 days 5 6 7 don't have either structure, so load posTrials in those...
        load(fullfile(basepath,['posTrials.mat']), 'posTrials'); posOn = posTrials{2}; posOff = posTrials{1};
    end
    q = sortrows([posOn; posOff]); behavior.position.x = q(:,2)';
    behavior.position.y = ones(size(q(:,2)))'; behavior.positionTrials{1} = posOn; behavior.positionTrials{2} = posOff;
    isBoundary = zscore(diff(posOn(:,1)))>1; trialsOn = [posOn([true; isBoundary],1) posOn([isBoundary;true],1)];
    isBoundary = zscore(diff(posOff(:,1)))>1; trialsOff = [posOff([true; isBoundary],1) posOff([isBoundary;true],1)];
    trials = Group(trialsOn,trialsOff); trials(:,3) = 2-trials(:,3); behavior.trials = trials(:,1:2); behavior.trialID = trials(:,[3 3]);
end
try run = behavior.run;
catch
    ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:)); t = behavior.timestamps(ok);
    speed = LinearVelocity([behavior.timestamps(ok)' behavior.position.x(ok)' behavior.position.y(ok)'],5);
    allPositions = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}]);
    [up,down] = ZeroCrossings([allPositions(:,1),allPositions(:,2)-0.5]);
    midpoints = allPositions(up|down,1);
    topSpeed = interp1(speed(:,1),speed(:,2),midpoints);
    threshold = median(topSpeed)/10; % the threshold is 10% of the peak speed
    run = t(FindInterval(speed(:,2)>threshold));
    run = ConsolidateIntervals(run,'epsilon',0.01);
    [in,w] = InIntervals(behavior.timestamps(:),run);
    peak = Accumulate(w(in),behavior.speed(in)','mode','max');
    % remove outliers (data in between sessions gives outlier speeds)
    [~,isOutlier] = RemoveOutliers(peak);
    % remove run epochs that don't reach the speed threshold
    run(peak<0.1 | isOutlier,:) = [];
    run(IntervalsIntersect(run,sleep),:) = [];
    % Make sure our each run epoch is confined within the same trial: no run epoch should begin in one trial and continue over the next trial
    run = SubtractIntervals(run,bsxfun(@plus,sort(behavior.trials(:)),[-1 1]*0.00001));
    run(diff(run,[],2)<0.5,:) = []; % remove run epochs lasting for less than 0.5s
end
onTrials = behavior.trials(behavior.trialID(:,2)==1,:); % Can promises that behavior.trialID(:,2)==1 is always stim and 0 is always nonstim trials
offTrials = behavior.trials(behavior.trialID(:,2)==0,:);
stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
stimIntervals = SubtractIntervals(stimIntervals,SubtractIntervals([0 Inf],onTrials)); % make sure stimIntervals only take place during onTrials
nonstimIntervals = SubtractIntervals(run,stim);
nonstimIntervals = SubtractIntervals(nonstimIntervals,SubtractIntervals([0 Inf],offTrials)); % make sure nonstimIntervals only take place during offTrials

id = spikes(:,2);
try load(fullfile(basepath,'thetaCyclesTask.mat'),'cycles');
catch load(fullfile(basepath,'thetaCyclesTrack.mat'),'cycles');
end

periods = {stimIntervals; nonstimIntervals};

%%
try
    for k=1:2,
        tic;
        period = periods{k};
        theseCycles = Restrict([cycles(:,1) cycles(:,2)],period);
        cycletime = RelativeTimes(spikes(:,1),theseCycles,[0 1]);
        [~,w] = InIntervals(spikes(:,1),theseCycles);
        cycletime = cycletime+w-1;
        ok = ~isnan(cycletime);
        for i=1:length(spikesCell),
            %         [spectrogram,t,f] = PointSpectrogram(cycletime(id==i & ok),'range',[0 3],'step',1,'window',80,'pad',2);
            try
            [spectrogram,t,f] = PointSpectrogram(cycletime(id==i & ok),'range',[0 3],'step',1,'window',floor(max(cycletime(id==i & ok)))-1,'pad',0);
            SpectrogramCell{1,k}(i,:) = nanmean(spectrogram',1);
            catch 
                if i>1,SpectrogramCell{1,k}(i,:) = nan; end
                display(['caught ' num2str(i)]);
            end
            if k==1, disp(['Stim: done with stim unit ' num2str(i) '/' num2str(length(spikesCell)) '; estimated time remaining=' num2str(toc/i*(length(spikesCell)-i))]); end
            if k==2, disp(['Nonstim: done with unit ' num2str(i) '/' num2str(length(spikesCell)) '; estimated time remaining=' num2str(toc/i*(length(spikesCell)-i))]); end
        end
    end
    specStim = SpectrogramCell{1}; specNonStim = SpectrogramCell{2};
catch
    keyboard
end
try
    save(['M:\home\raly\Documents\code\new\tempFiles\' name],'specStim','specNonStim','f');
    display('saved!');
end

% transform to a common number of bins
raw = specStim; f0 = f;
f = linspace(0,1,501)';
specStim = interp1(f0,raw,f);
raw = specNonStim;
specNonStim = interp1(f0,raw,f);





 