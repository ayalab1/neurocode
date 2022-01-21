function behavior = general_behavior_file(varargin)
% converts multiple tracking data types to the standard described in cellexplorer
% https://cellexplorer.org/datastructure/data-structure-and-format/#behavior
%
% This was writed to standardize xy coordinates and trials in several older data sets
% This is not the cleanest code, but takes care of most cases
%
% check extract_tracking below to preview methods. Can be further
% customized.
%
% Ryan H 2021

p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'fs',39.0625); % behavioral tracking sample rate (will detect fs for newer datasets)
addParameter(p,'force_overwrite',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'force_run',true); % run even if animal.behavior already exists
addParameter(p,'save_mat',true); % save animal.behavior.mat

parse(p,varargin{:});
basepaths = p.Results.basepath;
fs = p.Results.fs;
force_overwrite = p.Results.force_overwrite;
force_run = p.Results.force_run;
save_mat = p.Results.save_mat;

if ~iscell(basepaths)
    basepaths = {basepaths};
end

for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file') &&...
            ~force_run
        continue
    end
    disp(basepath)
    behavior = main(basepath,basename,fs,save_mat,force_overwrite);
end
end

function behavior = main(basepath,basename,fs,save_mat,force_overwrite)

if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file') &&...
        ~force_overwrite
    load([basepath,filesep,[basename,'.animal.behavior.mat']]);
end

[t,x,y,z,v,trials,units,source,linearized,fs,notes] =...
    extract_tracking(basepath,basename,fs);

load([basepath,filesep,[basename,'.session.mat']]);

behavior.sr = fs;
behavior.time = t';
behavior.position.x = x';
behavior.position.y = y';
behavior.position.z = z';
behavior.position.linearized = linearized';
behavior.position.units = units;
behavior.speed = v';
behavior.acceleration = [0,diff(behavior.speed)];
behavior.trials = trials;
behavior.notes = notes;
behavior.epochs = session.epochs;
behavior.processinginfo.date = date;
behavior.processinginfo.function = 'general_behavioral_file.mat';
behavior.processinginfo.source = source;

if save_mat
    save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
end
end

function [t,x,y,z,v,trials,units,source,linearized,fs,notes] =...
    extract_tracking(basepath,basename,fs)

t = [];
x = [];
y = [];
z = [];
v = [];
trials = [];
units = [];
source = [];
linearized = [];
notes = [];

% below are many methods on locating tracking data from many formats

% standard whl file xyxy format
if exist([basepath,filesep,[basename,'.whl']],'file')
    positions = load([basepath,filesep,[basename,'.whl']]);
    t = (0:length(positions)-1)'/fs;
    positions(positions == -1) = NaN;
    % find led with best tracking
    [x,y] = find_best_columns(positions,fs);
    units = 'cm';
    source = '.whl';
    if exist(fullfile(basepath,'trials.mat'),'file')
        trials = load(fullfile(basepath,'trials.mat'));
        trials = trials.trials;
    end
    if exist(fullfile(basepath,[basename,'-TrackRunTimes.mat']),'file')
        trials = load(fullfile(basepath,[basename,'-TrackRunTimes.mat']));
        trials = trials.trackruntimes;
    end
    
    % sometimes whl files are within fmat folder and have different name
elseif exist(fullfile(basepath,'fmat',[animalFromBasepath(basepath),basename,'.whl']),'file')
    
    positions = load(fullfile(basepath,'fmat',[animalFromBasepath(basepath),basename,'.whl']));
    t = (0:length(positions)-1)'/fs;
    positions(positions == -1) = NaN;
    try
        [x,y] = find_best_columns(positions,fs);
    catch
        x = positions(:,1);
        y = positions(:,2);
    end
    units = 'cm';
    source = '.whl';
    if exist(fullfile(basepath,'fmat','trials.mat'),'file')
        trials = load(fullfile(basepath,'fmat','trials.mat'));
        try
            trials = trials.trials;
        catch
            temp_trials = [];
            for name = fieldnames(trials)'
                temp_trials = [temp_trials;trials.(name{1})];
            end
            [~,idx] = sort(temp_trials(:,1));
            trials = temp_trials(idx,:);
        end
    end
    
    % sometimes whl files are within fmat folder and have different name
elseif exist(fullfile(basepath,'fmat',[basename,'.whl']),'file')
    positions = load(fullfile(basepath,'fmat',[basename,'.whl']));
    t = (0:length(positions)-1)'/fs;
    positions(positions == -1) = NaN;
    % find led with best tracking
    [x,y] = find_best_columns(positions,fs);
    units = 'cm';
    source = '.whl';
    if exist(fullfile(basepath,'fmat','trials.mat'),'file')
        trials = load(fullfile(basepath,'fmat','trials.mat'));
        trials = trials.trials;
    end
    
elseif ~isempty(dir(fullfile(basepath,'fmat', '*.whl')))
    filelist = dir(fullfile(basepath,'fmat', '*.whl'));
    positions = load(fullfile(filelist(1).folder,filelist(1).name));
    t = (0:length(positions)-1)'/fs;
    positions(positions == -1) = NaN;
    % find led with best tracking
    [x,y] = find_best_columns(positions,fs);
    units = 'cm';
    source = '.whl';
    if exist(fullfile(basepath,'fmat','trials.mat'),'file')
        trials = load(fullfile(basepath,'fmat','trials.mat'));
        trials = trials.trials;
    elseif exist(fullfile(basepath,[basename,'.trials.mat']),'file')
        trials = load(fullfile(basepath,[basename,'.trials.mat']));
        temp_trials = [];
        for name = fieldnames(trials)'
            temp_trials = [temp_trials;trials.(name{1})];
        end
        [~,idx] = sort(temp_trials(:,1));
        trials = temp_trials(idx,:);
    end
    % postTrials format, processed linearized data
elseif exist([basepath,filesep,['posTrials.mat']],'file')
    load([basepath,filesep,['posTrials.mat']]);
    positions = [posTrials{1};posTrials{2}];
    [~,idx] = sort(positions(:,1));
    positions = positions(idx,:);
    t = positions(:,1);
    x = [];
    y = [];
    linearized = positions(:,2);
    units = 'normalize';
    source = 'posTrials.mat';
    fs = 1/mode(diff(t));
    if exist(fullfile(basepath,[basename,'.trials.mat']),'file')
        trials = load(fullfile(basepath,[basename,'.trials.mat']));
        temp_trials = [];
        for name = fieldnames(trials)'
            temp_trials = [temp_trials;trials.(name{1})];
        end
        [~,idx] = sort(temp_trials(:,1));
        trials = temp_trials(idx,:);
    end
    
    % posTrials is sometimes moved
elseif exist([basepath,filesep,['oldfiles\posTrials.mat']],'file')
    load([basepath,filesep,['oldfiles\posTrials.mat']]);
    positions = [posTrials{1};posTrials{2}];
    [~,idx] = sort(positions(:,1));
    positions = positions(idx,:);
    t = positions(:,1);
    x = [];
    y = [];
    linearized = positions(:,2);
    units = 'normalize';
    source = 'posTrials.mat';
    fs = 1/mode(diff(t));
    
    if exist(fullfile(basepath,[basename,'.trials.mat']),'file')
        trials = load(fullfile(basepath,[basename,'.trials.mat']));
        temp_trials = [];
        for name = fieldnames(trials)'
            temp_trials = [temp_trials;trials.(name{1})];
        end
        [~,idx] = sort(temp_trials(:,1));
        trials = temp_trials(idx,:);
    end
        
    % .position.behavior file with x,y,linear and more
elseif exist([basepath,filesep,[basename,'.position.behavior.mat']],'file')
    load([basepath,filesep,[basename,'.position.behavior.mat']])
    t = position.timestamps;
    x = position.position.x;
    y = position.position.y;
    linearized = position.position.lin;
    units = position.units;
    
    if position.units == "m"
        x = x*100;
        y = y*100;
        linearized = linearized*100;
        units = 'cm';
    end
    source = 'position.behavior.mat';
    fs = 1/mode(diff(t));
    
    if exist([basepath,filesep,['position_info.mat']],'file')
        load([basepath,filesep,['position_info.mat']])
         trials = [cellfun(@(x) min(x),pos_inf.ts_ep),...
            cellfun(@(x) max(x),pos_inf.ts_ep)];
    end
    
    % position_info files have xy and linearized data
elseif exist([basepath,filesep,['position_info.mat']],'file')
    load([basepath,filesep,['position_info.mat']])
    t = pos_inf.ts';
    x = pos_inf.x;
    y = pos_inf.y;
    linearized = pos_inf.lin_pos;
    trials = [cellfun(@(x) min(x),pos_inf.ts_ep),...
        cellfun(@(x) max(x),pos_inf.ts_ep)];
    units = 'cm';
    source = 'position_info.mat';
    fs = 1/mode(diff(t));
    
    %  _TXVt files have time, x, v, and trials
elseif exist([basepath,filesep,[basename,'_TXVt.mat']],'file')
    load([basepath,filesep,[basename,'_TXVt.mat']])
    t = TXVt(:,1);
    linearized = TXVt(:,2);
    x = TXVt(:,2);
    y = [];
    for trial_n = unique(TXVt(:,4))'
        trial_ts = TXVt(TXVt(:,4) == trial_n,1);
        trials = [trials;[min(trial_ts),max(trial_ts)]];
    end
    units = 'cm';
    source = '_TXVt.mat';
    fs = 1/mode(diff(t));
    
elseif exist([basepath,filesep,[basename,'.tracking.behavior.mat']],'file')
    load([basepath,filesep,[basename,'.tracking.behavior.mat']])
    t = tracking.timestamps;
    x = tracking.position.x * 100;
    y = tracking.position.z * 100;
    z = tracking.position.y * 100;
    notes = "z to y and y to z";
    units = 'cm';
    source = '.tracking.behavior.mat';
    
    fs = 1/mode(diff(t));
    
    if ~isempty(dir([basepath,filesep,[basename,'.*Trials.mat']]))
        filelist = dir([basepath,filesep,[basename,'.*Trials.mat']]);
        for file = filelist'
            load(fullfile(file.folder,file.name))
            if ~isempty(trials)
                trials = trials.int;
                break
            end
        end
    end
    % sometimes whl files have yet to be concatenated and are in subfolders
    % elseif ~isempty(dir(fullfile(basepath, '**\*.whl')))
    %     filelist = dir(fullfile(basepath, '**\*.whl'));
    %     for file = filelist'
    %         disp(file)
    %     end
    %
else
    warning('write another loader here')
    return
end

% trials can sometimes have extra columns
if size(trials,2) > 2
    trials = trials(:,1:2);
end
% to help find if trials are index instead of sec
isaninteger = @(x)isfinite(x) & x==floor(x);
% check if trials are integers, if so, they are index instead of sec
if all(isaninteger(trials)) & ~isempty(trials)
    trials = t(trials);
end
% if the max trial is greater than the available time, they are index
% if max(trials(:))-max(t) > 10 & all(isaninteger(trials)) & ~isempty(trials)
%     trials = t(trials);
% end

% get velocity
try
    try
        v = LinearVelocity([t,x,y]);
        v = v(:,2);
    catch
        [~,idx,idx1] = unique(t);
        v = LinearVelocity([t(idx),x(idx),y(idx)]);
        v = v(idx1,:);
        v = v(:,2);
    end
catch
    v = LinearVelocity([t,linearized,linearized*0]);
    v = v(:,2);
end

end

% find led with best tracking
function [x,y] = find_best_columns(positions,fs)
for col = 1:size(positions,2)
    x_test = medfilt1(positions(:,col),round(fs/2),'omitnan');
    R(col) = corr(positions(:,col),x_test, 'rows','complete');
end
[~,idx] = max([mean(R(1:2)), mean(R(3:4))]);
columns{1} = [1,2];
columns{2} = [3,4];
x = positions(:,columns{idx}(1));
y = positions(:,columns{idx}(2));
%     [~,idx] = min([sum(isnan(positions(:,1))),sum(isnan(positions(:,3)))]);

end