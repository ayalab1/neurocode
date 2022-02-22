function tracking = process_and_sync_dlc(varargin)
% Unpacks DLC CSV within subfolders
%
% Run this after you have exported deeplabcut csv results
%
% TODO: 
%       -make multi tracking points available (multi-points on animal)
%       -make option to just load and format csv without sync
%       -make more robust to file exist assumptions 
%           (ex. MergePoints.events.mat, digitalIn.events.mat)
%       -make option for dlc csv in basepath 
%       -make compatible with multi-animal tracking
%
%
% Ryan H 2022

p = inputParser;
p.addParameter('basepath',pwd,@isfolder); 
p.addParameter('primary_coords',1,@isnumeric); % which tracker point do you want
p.addParameter('likelihood',.95,@isnumeric); % tracking quality thres [0-1]
p.addParameter('pulses_delta_range',0.01,@isnumeric); % range for ttls

p.parse(varargin{:});
primary_coords = p.Results.primary_coords;
basepath = p.Results.basepath;
likelihood = p.Results.likelihood;
pulses_delta_range = p.Results.pulses_delta_range;

basename = basenameFromBasepath(basepath);
            
if exist(fullfile(basepath,[basename,'.MergePoints.events.mat']),'file')
    load(fullfile(basepath,[basename,'.MergePoints.events.mat']));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        if ~isempty(dir(fullfile(basepath,MergePoints.foldernames{ii},'*DLC*csv')))
            % locate file
            file = dir(fullfile(basepath,MergePoints.foldernames{ii},'*DLC*csv'));
            
            video_file = dir(fullfile(file.folder,'*.avi'));
            obj = VideoReader(fullfile(video_file.folder,video_file.name));
            fs = obj.FrameRate;
            
            % load csv with proper header
            opts = detectImportOptions(fullfile(file.folder,file.name),'NumHeaderLines',2);
            df = readtable(fullfile(file.folder,file.name),opts);
            % get names of fields, these will be as long as tracker points
            % used times 3 because [x,y,likelihood]
            field_names = fields(df);
            % locate columns with [x,y,likelihood]
            x_col = find(contains(field_names,'x'));
            y_col = find(contains(field_names,'y'));
            likelihood_col = find(contains(field_names,'likelihood'));
            % filter out bad tracker points by likelihood thres
            for i = 1:length(x_col)
                idx = df{:,likelihood_col(i)} < likelihood;
                df{idx,x_col(i)} = NaN;
                df{idx,y_col(i)} = NaN;
            end
            ts = df{:,1}/fs;
            x = df{:,x_col(primary_coords)};
            y = df{:,y_col(primary_coords)};
            
            tempTracking{count} = sync_ttl(file.folder,x,y,ts,fs,pulses_delta_range);
            trackFolder(count) = ii;
            count = count + 1;
        end
    end
end
    
% Concatenate and sync timestamps
ts = []; subSessions = []; maskSessions = [];
if exist(fullfile(basepath,[basename,'.MergePoints.events.mat']),'file')
    load(fullfile(basepath,[basename,'.MergePoints.events.mat']));
    for ii = 1:length(trackFolder)
        if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
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

% Concatenating tracking fields...
x = []; y = []; folder = []; samplingRate = []; description = [];
for ii = 1:size(tempTracking,2)
    x = [x; tempTracking{ii}.position.x];
    y = [y; tempTracking{ii}.position.y];
    folder{ii} = tempTracking{ii}.folder;
    samplingRate = [samplingRate; tempTracking{ii}.samplingRate];
    description{ii} = tempTracking{ii}.description;
end
    
tracking.position.x = x;
tracking.position.y = y;
tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = ts;
tracking.events.subSessions =  subSessions;
tracking.events.subSessionsMask = maskSessions;
end

function [tracking] = sync_ttl(folder,x,y,ts,fs,pulses_delta_range)

load(fullfile(folder,'digitalIn.events.mat'))
bazlerTtl = digitalIn.timestampsOn{1,1};

%check for extra pulses of much shorter distance than they should
extra_pulses = diff(bazlerTtl)<((1/fs)-(1/fs)*pulses_delta_range);
bazlerTtl(extra_pulses) = [];

basler_intan_diff = length(bazlerTtl) - length(x);

[x,y,ts,bazlerTtl] = match_basler_frames_to_ttl(bazlerTtl,basler_intan_diff,x,y,ts,fs);

[~,folder_name] = fileparts(folder);
tracking.position.x =x;
tracking.position.y = y;
tracking.timestamps = bazlerTtl;
tracking.originalTimestamps = ts;
tracking.folder = folder_name;
tracking.samplingRate = fs;
tracking.description = '';
end


function [x,y,t,bazlerTtl] = match_basler_frames_to_ttl(bazlerTtl,basler_intan_diff,x,y,t,fs)

% match basler frames con ttl pulses
if (length(bazlerTtl) == length(x)) || abs(basler_intan_diff)<=2 %assumes 1 frame could be cut at 0 and 1 frame at end
    disp('N of frames match!!');
elseif basler_intan_diff>0 && abs(basler_intan_diff)<fs
    disp([num2str(abs(length(bazlerTtl) - length(x))) ' of frames dont match, probably at the end of the recording']);
    bazlerTtl = bazlerTtl(1:length(x));
elseif basler_intan_diff<0 && abs(basler_intan_diff)<fs
    disp([num2str(abs(length(bazlerTtl) - length(x))) ' of frames dont match, probably at the beggining of the recording']);
    x = x(1:length(bazlerTtl));
    y = y(1:length(bazlerTtl));
elseif basler_intan_diff<0 && abs(basler_intan_diff)>fs
    disp([num2str(abs(length(x) - length(bazlerTtl)))...
        ' video frames without TTL... was the recording switched off before the camera? Cutting positions accordingly...']);
    x = x(1:length(bazlerTtl));
    y = y(1:length(bazlerTtl));
elseif abs(basler_intan_diff)>2*fs
    warning('More than 2 seconds missalignment in total in this session...will adjust to the closer one...');
    if basler_intan_diff>0
        bazlerTtl = bazlerTtl(1:length(x));
    else
        x = x(1:length(bazlerTtl));
        y = y(1:length(bazlerTtl));
    end
elseif isempty(bazlerTtl)
    bazlerTtl = t;
else
    warning('Unnoticed problem with Camera/Intan... I would go back and check both step by step');
end
end
