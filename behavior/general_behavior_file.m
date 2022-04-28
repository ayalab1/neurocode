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
% Currently compatible with the following sources:
%   .whl, posTrials.mat, basename.posTrials.mat, position.behavior.mat, position_info.mat,
%   _TXVt.mat, tracking.behavior.mat, Tracking.behavior.mat, DeepLabCut
%
% Ryan H 2021

p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'fs',39.0625); % behavioral tracking sample rate (will detect fs for newer datasets)
addParameter(p,'force_overwrite',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'force_run',true); % run even if animal.behavior already exists
addParameter(p,'save_mat',true); % save animal.behavior.mat
addParameter(p,'primary_coords_dlc',1); % deeplabcut tracking point to extract (extracts all, but main x and y will be this)
addParameter(p,'likelihood_dlc',.95); % deeplabcut likelihood threshold

parse(p,varargin{:});
basepaths = p.Results.basepath;
fs = p.Results.fs;
force_overwrite = p.Results.force_overwrite;
force_run = p.Results.force_run;
save_mat = p.Results.save_mat;
primary_coords_dlc = p.Results.primary_coords_dlc;
likelihood_dlc = p.Results.likelihood_dlc;

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
    behavior = main(basepath,basename,fs,save_mat,force_overwrite,...
        primary_coords_dlc,likelihood_dlc);
end
end

function behavior = main(basepath,basename,fs,save_mat,force_overwrite,primary_coords,likelihood)

if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file') &&...
        ~force_overwrite
    load([basepath,filesep,[basename,'.animal.behavior.mat']]);
end

[t,x,y,z,v,trials,units,source,linearized,fs,notes,extra_points,stateNames,states] =...
    extract_tracking(basepath,basename,fs,primary_coords,likelihood);

load([basepath,filesep,[basename,'.session.mat']]);

behavior.sr = fs;
behavior.timestamps = t';
behavior.position.x = x';
behavior.position.y = y';
behavior.position.z = z';
behavior.position.linearized = linearized';
behavior.position.units = units;
behavior.speed = v';
behavior.acceleration = [0,diff(behavior.speed)];
behavior.trials = trials;
behavior.states = states;
behavior.stateNames = stateNames;
behavior.notes = notes;
behavior.epochs = session.epochs;
behavior.processinginfo.date = date;
behavior.processinginfo.function = 'general_behavioral_file.mat';
behavior.processinginfo.source = source;

if ~isempty(extra_points)
    for field = fieldnames(extra_points)'
        field = field{1};
        behavior.position.(field) = extra_points.(field)';
    end
end

if save_mat
    save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
end
end

function [t,x,y, z,v,trials,units,source,linearized,fs,notes,extra_points,...
    stateNames,states] = extract_tracking(basepath,basename,fs,primary_coords,likelihood)

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
extra_points = [];
stateNames = [];
states = [];
% below are many methods on locating tracking data from many formats


% search for DLC csv within basepath and subdirs, but not kilosort folder
%       (takes too long)
if exist(fullfile(basepath,[basename,'.MergePoints.events.mat']),'file')
    
    load(fullfile(basepath,[basename,'.MergePoints.events.mat']))
    for k = 1:length(MergePoints.foldernames)
        dlc_flag(k) = ~isempty(dir(fullfile(basepath,MergePoints.foldernames{k},'*DLC*.csv')));
    end
    files = dir(basepath);
    files = files(~contains({files.name},'Kilosort'),:);
    dlc_flag(k+1) = ~isempty(dir(fullfile(files(1).folder,'*DLC*.csv')));
else
    dlc_flag = false;
end
if any(dlc_flag)
    disp('detected deeplabcut')
    [tracking,field_names] = process_and_sync_dlc('basepath',basepath,...
        'primary_coords',primary_coords,...
        'likelihood',likelihood);
    
    t = tracking.timestamps;
    fs = 1/mode(diff(t));
    
    x = tracking.position.x(:,primary_coords);
    y = tracking.position.y(:,primary_coords);
    
    % multiple tracking points will likely exist, extract here
    x_col = field_names(contains(field_names,'x'));
    y_col = field_names(contains(field_names,'y'));
    extra_points = struct();
    for i = 1:length(x_col)
        extra_points.([x_col{i},'_point']) = tracking.position.x(:,i);
        extra_points.([y_col{i},'_point']) = tracking.position.y(:,i);
    end
    
    units = 'pixels';
    source = 'deeplabcut';
    
    if length(t) > length(x)
        t = t(1:length(x));
    elseif length(x) > length(t)
        x = x(1:length(t));
        y = y(1:length(t));
        % adjust other tracker points
        for name = fields(extra_points)'
            extra_points.(name{1}) = extra_points.(name{1})(1:length(t));
        end
    end
    
    if isfield(tracking, 'events')
        if isfield(tracking.events,'subSessions')
            trials = tracking.events.subSessions;
        end
    end
    notes = ['primary_coords: ',num2str(primary_coords),...
        ', likelihood: ',num2str(likelihood)];
    notes = {notes,tracking.notes};
    
    % standard whl file xyxy format
elseif exist([basepath,filesep,[basename,'.whl']],'file')
    disp('detected .whl')
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
    disp('detected .whl')
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
    disp('detected .whl')
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
    disp('detected .whl')
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
    disp('detected posTrials.mat')
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
elseif exist([basepath,filesep,[basename,'.posTrials.mat']],'file')
    disp('detected basename.posTrials.mat')
    load([basepath,filesep,[basename,'.posTrials.mat']]);
    
    positions = [];
    linearized = [];
    states_temp = [];
    trials = [];
    for ep = 1:length(posTrials)
        positions = [positions;posTrials{ep}.pos];
        linearized = [linearized;posTrials{ep}.linpos];
        states_temp = [states_temp;repmat({posTrials{ep}.type},length(posTrials{ep}.linpos),1)];
        trials = [trials;posTrials{ep}.int];
    end
    
    [~,idx] = sort(positions(:,1));
    positions = positions(idx,:);
    
    [~,idx] = sort(linearized(:,1));
    linearized = linearized(idx,1);
    
    states_temp = states_temp(idx,:);
    
    stateNames = unique(states_temp)';
    states = zeros(1,length(states_temp));
    for i = 1:length(stateNames)
        states(contains(states_temp,stateNames{i})) = i;
    end
    
    t = positions(:,1);
    x = positions(:,2);
    y = positions(:,3);
    
    units = 'normalize';
    source = 'basename.posTrials.mat';
    
    % posTrials is sometimes moved
elseif exist([basepath,filesep,['oldfiles\posTrials.mat']],'file')
    disp('detected posTrials.mat')
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
    disp('detected position.behavior.mat')
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
    disp('detected position_info.mat')
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
    disp('detected _TXVt.mat')
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
    disp('detected tracking.behavior.mat')
    load([basepath,filesep,[basename,'.tracking.behavior.mat']])
    t = tracking.timestamps;
    fs = 1/mode(diff(t));
    
    if isfield(tracking.position,'x') && isfield(tracking.position,'y') && isfield(tracking.position,'z')
        x = tracking.position.x * 100;
        y = tracking.position.z * 100;
        z = tracking.position.y * 100;
        notes = "z to y and y to z";
        units = 'cm';
        source = '.tracking.behavior.mat';
        
    elseif isfield(tracking.position,'x1') && isfield(tracking.position,'y1')
        positions = [tracking.position.x1,tracking.position.y1,tracking.position.x2,tracking.position.y2];
        [x,y] = find_best_columns(positions,fs);
        if range(x) <= 1
            units = 'normalized';
        elseif range(x) > 1
            units = 'pixels';
        end
        source = '.Tracking.behavior.mat';
        if length(t) > length(x)
            t = t(1:length(x));
            warning('Different number of ts and coords! Check data source to verify')
        end
    end
    if isfield(tracking, 'events')
        if isfield(tracking.events,'subSessions')
            trials = tracking.events.subSessions;
        end
    end
    
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
    warning('No video detected...')
    disp('attempting to add ttls from digitalIn')
    if exist(fullfile(basepath,[basename,'.MergePoints.events.mat']),'file')
        load(fullfile(basepath,[basename,'.MergePoints.events.mat']));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            tempTracking{count} = sync_ttl(basepath,MergePoints.foldernames{ii});
            trackFolder(count) = ii;
            count = count + 1;
        end
    end
    % Concatenate and sync timestamps
    ts = []; subSessions = []; maskSessions = [];
    if exist(fullfile(basepath,[basename,'.MergePoints.events.mat']),'file')
        load(fullfile(basepath,[basename,'.MergePoints.events.mat']));
        for ii = 1:length(trackFolder)
            if strcmpi(fullfile(basepath,MergePoints.foldernames{trackFolder(ii)}),tempTracking{ii}.folder)
                sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
                subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
                maskSessions = [maskSessions; ones(size(sumTs))*ii];
                ts = [ts; sumTs];
            else
                error('Folders name does not match!!');
            end
        end
    else
        warning('No MergePoints file found. Concatenating timestamps...');
        for ii = 1:length(trackFolder)
            sumTs = max(ts)+ tempTracking{ii}.timestamps;
            subSessions = [subSessions; [sumTs(1) sumTs(end)]];
            ts = [ts; sumTs];
        end
    end
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

function tracking = sync_ttl(basepath,folder)

if ~exist(fullfile(basepath,folder,'digitalIn.events.mat'),'file')
    digitalIn = getDigitalIn('all','folder',fullfile(basepath,folder));
end
load(fullfile(fullfile(basepath,folder),'digitalIn.events.mat'))

Len = cellfun(@length, digitalIn.timestampsOn, 'UniformOutput', false);
[~,idx] = max(cell2mat(Len));
bazlerTtl = digitalIn.timestampsOn{idx};
fs = mode(diff(bazlerTtl));
%check for extra pulses of much shorter distance than they should
extra_pulses = diff(bazlerTtl)<((1/fs)-(1/fs)*0.01);
bazlerTtl(extra_pulses) = [];

[~,folder_name] = fileparts(folder);
tracking.timestamps = bazlerTtl;
tracking.folder = fullfile(basepath,folder);
tracking.samplingRate = fs;
tracking.description = '';
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

function [digitalIn] = getDigitalIn(ch,varargin)
% [pul, val, dur] = getDigitalIn(d,varargin)
%
% Find digital In pulses
%
% INPUTS
% ch            Default all.
% <OPTIONALS>
% fs            Sampling frequency (in Hz), default 30000, or try to
%               recover for rhd
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)
% filename      File to get pulses from. Default, digitalin.dat file with folder
%               name in current directory
%
%
% OUTPUTS
%               digitalIn - events struct with the following fields
% ints          C x 2  matrix with pulse times in seconds. First column of C
%               are the beggining of the pulses, second column of C are the end of
%               the pulses.
% dur           Duration of the pulses. Note that default fs is 30000.
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag
%
% MV-BuzsakiLab 2019
% Based on Process_IntanDigitalChannels by P Petersen

% Parse options
if exist('ch') ~= 1
    ch = 'all';
end

p = inputParser;
addParameter(p,'fs',[],@isnumeric)
addParameter(p,'offset',0,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',5,@isnumeric)
addParameter(p,'folder',pwd,@isfolder)

parse(p, varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
folder = p.Results.folder;

if ~isempty(dir(fullfile(folder,'*.xml')))
    %sess = bz_getSessionInfo(pwd,'noPrompts',true);
    sess = getSession('basepath',fileparts(folder));
end
if ~isempty(dir(fullfile(folder,'*DigitalIn.events.mat')))
    disp('Pulses already detected! Loading file.');
    file = dir(fullfile(folder,'*DigitalIn.events.mat'));
    load(fullfile(file.folder,file.name));
    return
end

if isempty(filename)
    filename=dir(fullfile(folder,'digitalIn.dat'));
    filename = filename.name;
elseif exist('filename','var')
    disp(['Using input: ',filename])
else
    disp('No digitalIn file indicated...');
end

try [amplifier_channels, notes, aux_input_channels, spike_triggers,...
        board_dig_in_channels, supply_voltage_channels, frequency_parameters,board_adc_channels] =...
        read_Intan_RHD2000_file_bz('basepath',folder);
    fs = frequency_parameters.board_dig_in_sample_rate;
catch
    disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
end

disp('Loading digital channels...');
m = memmapfile(fullfile(folder,filename),'Format','uint16','writable',false);
digital_word2 = double(m.Data);
clear m
Nchan = 16;
Nchan2 = 17;
for k = 1:Nchan
    tester(:,Nchan2-k) = (digital_word2 - 2^(Nchan-k))>=0;
    digital_word2 = digital_word2 - tester(:,Nchan2-k)*2^(Nchan-k);
    test = tester(:,Nchan2-k) == 1;
    test2 = diff(test);
    pulses{Nchan2-k} = find(test2 == 1);
    pulses2{Nchan2-k} = find(test2 == -1);
    data(k,:) = test;
end
digital_on = pulses;
digital_off = pulses2;
disp('Done!');

for ii = 1:size(digital_on,2)
    if ~isempty(digital_on{ii})
        % take timestamp in seconds
        digitalIn.timestampsOn{ii} = digital_on{ii}/fs;
        digitalIn.timestampsOff{ii} = digital_off{ii}/fs;
        
        % intervals
        d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
        d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
        d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        digitalIn.ints{ii} = d;
        digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % durantion
        
        clear intsPeriods
        intsPeriods(1,1) = d(1,1); % find stimulation intervals
        intPeaks =find(diff(d(1,:))>lag);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = d(2,intPeaks(jj));
            intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
        end
        intsPeriods(end,2) = d(2,end);
        digitalIn.intsPeriods{ii} = intsPeriods;
    end
end

if exist('digitalIn','var')==1
    xt = linspace(0,size(data,2)/fs,size(data,2));
    data = flip(data);
    data = data(1:size(digitalIn.intsPeriods,2),:);
    
    h=figure;
    imagesc(xt,1:size(data,2),data);
    xlabel('s'); ylabel('Channels'); colormap gray
    mkdir(fullfile(folder,'Pulses'));
    saveas(h,fullfile(folder,'Pulses','digitalIn.png'));
    
    save(fullfile(folder,'digitalIn.events.mat'),'digitalIn');
else
    digitalIn = [];
end

end