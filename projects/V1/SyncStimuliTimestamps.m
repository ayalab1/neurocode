
folder = pwd;

[parentFolder,dayName] = fileparts(folder);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

% Load saved data
q = dir('task*.times');
matrix = [];
for i=1:length(q)
    filename = q(i).name;
    data = dlmread(filename);
    % convert data to the time of the merged datfile:
    ok = ~isnan(data(:,1));
    if size(data,2)>4, continue; end % ignore Ymaze file
    data(:,end+1) = i;
    matrix = [matrix; data];
end
matrix(isnan(matrix(:,1)),:) = [];

try
    for i=1:length(session.epochs)
        isVS(i,1) = any(strfind(session.epochs{i}.name,'isual'));
    end
    if any(isVS)
        VS = session.epochs(find(isVS)); for i=1:length(VS), VS{i} = [VS{i}.startTime  VS{i}.stopTime]; end
        VS = cell2mat(VS');
    else error('');
    end
catch
    VS = [0 Inf]; % no restriction
end

% Use digital input if any
try
    %     exist(fullfile(parentFolder,dayName,'digitalIn.events.mat'),'file')
    load(fullfile(parentFolder,dayName,'digitalIn.events.mat'),'digitalIn');
    on = digitalIn.timestampsOn{end};
    off = digitalIn.timestampsOff{end};
    if off(1)<on(1), off(1) = []; end
    if off(end)<on(end), off(end+1) = nan; end
    timesDigital = [on off];
    timesDigital = reshape(timesDigital',1,[])'; % vectorize
    % perhaps do: timesDigital = Restrict(timesDigital,VS); % visual stimuli to avoid the Y maze stuff

    % convert matrix time to seconds:
    t = (matrix(:,1)-matrix(1))*24*3600;
    on = t(strfind(matrix(:,2)'==-2,[1 0])+1,1);
    off = t(strfind(matrix(:,2)'>-2,[1 0])+1,1);
    if off(1)<on(1), off(1) = []; end
    if off(end)<on(end), off(end+1) = nan; end
    timesMatrix = [on off];
    timesMatrix = reshape(timesMatrix',1,[])'; % vectorize
    ok = ~isnan(timesDigital) & ~isnan(timesMatrix);

    newTimes = interp1(timesMatrix(ok),timesDigital(ok),matrix(:,1));
    % extrapolate missing timestamps
    missing = find(isnan(newTimes)); ok = find(~isnan(newTimes)); closest = ok(FindClosest(ok,missing));
    newTimes(missing) = newTimes(closest) + (matrix(missing,1)-(matrix(closest,1)));
    matrix(:,1) = newTimes;
catch
    display('Mismatch between matlab matrix and digital inputs. Ignoring digital inputs...');
    load(fullfile(folder,[dayName '.MergePoints.events.mat']),'MergePoints');
    for i=1:length(MergePoints.foldernames)
        computerTime(i,1) = dir(fullfile(folder,MergePoints.foldernames{i},'amplifier.dat')).datenum;
    end
    for i=1:max(matrix(:,end))
        indices = find(matrix(:,end)==i);
        index = FindClosest(computerTime,matrix(indices(end),1),'higher'); % which subsession were the stimuli presented
        timeToEnd = (computerTime(index)-matrix(indices,1))*24*3600; % convert days to seconds
        saved{i,1} = timeToEnd; index
        datTime = MergePoints.timestamps(index,2)-timeToEnd;
        matrix(indices,1) = datTime;
    end
end

dlmwrite('stimuli_synced.times',matrix,'precision','%.6f');


