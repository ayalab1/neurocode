function [spikeT] = importSpikes(varargin)
% [spikes] = importSpikes(varargin)
% Import into the workspace spiketimes for selected units taking as input a
% pre-computed buzcode spikes.cellinfo structure 

% INPUTS
%    basepath        -path to recording where spikes.cellinfo is. Default: pwd
%    UID             -vector subset of UID's to load. Default: all.
%    channel         -specific channels (neuroscope indexing) to use. Default: all

% These inputs require CellExplorer cell_metrics pre calculated. NOT IMPLEMENTED YET
%    brainRegion     -string region ID to load neurons from specific region
%    cellType        -cell type to load
%    sleepState      -string sleep state to keep spikes falling within the
%                       specified sleep state. 
% *** If inputting multiple regions, types, or states, ensure the input 
%     follows the format: ["State1","State2"] ***

% OUTPUT
%    spikeT - cellinfo struct with the following fields
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .sessionName             

% TODO:
% Do we want more output fields?

% AntonioFR, 8/20; Lindsay K, 11/21

%% Parse inputs 

p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'brainRegion','',@isstring); 
addParameter(p,'cellType','',@isstring); 
addParameter(p,'sleepState','',@isstring);
addParameter(p,'UID',[],@isvector);
addParameter(p,'spikes',[],@isstruct);  
addParameter(p,'session',[],@isstruct);  
addParameter(p,'channel',[],@isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
region = p.Results.brainRegion;
type = p.Results.cellType;
state = p.Results.sleepState;
UID = p.Results.UID;
spikes = p.Results.spikes;
session = p.Results.session;
channel = p.Results.channel;

%% Load spikes 
basename = basenameFromBasepath(basepath);

if isempty(spikes) && exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'file') 
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
end

spikeT.UID = [];
spikeT.times = {};
            
%% Remove bad channels before we start
load([basepath, filesep, basename, '.cell_metrics.cellinfo.mat']);
if isfield(cell_metrics, 'tags')
    if isfield(cell_metrics.tags, 'Bad')
    ct = 1;
    for i = 1:length(spikes.UID)
        if isempty(find(cell_metrics.tags.Bad==spikes.UID(i),1))
            spikeT.UID(ct) = spikes.UID(i);
            spikeT.times{ct} = spikes.times{i};
            ct = ct+1;
        end
    end
    clear ct
    else
    spikeT.UID = spikes.UID;
    spikeT.times = spikes.times;        
    end
else
spikeT.UID = spikes.UID;
spikeT.times = spikes.times;
end

%% Output structure
if ~isempty(UID)
   spikeT.UID = spikes.UID(UID);
   spikeT.times = spikes.times(UID);
elseif ~isempty(channel)
    setUn = [];
    for i = 1:length(channel)
        setUn = [setUn find(spikes.maxWaveformCh == channel(i))];
    end
    setUn = sort(setUn);
    spikeT.UID = spikes.UID(setUn);
    spikeT.times = spikes.times(setUn);
end

% only one region will be assigned to each cell
if ~isempty(region)
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
    keepUID = []; keepTimes = [];
    for i = 1:length(region)
        tempUID = []; tempTimes = [];
        tempUID = spikes.UID(contains(cell_metrics.brainRegion, region(i)));
        tempTimes = spikes.times(contains(cell_metrics.brainRegion, region(i)));
        [~,useIndTemp,useInd] = intersect(tempUID, spikeT.UID);
        keepUID = [keepUID spikeT.UID(useInd)];
        for j = 1:length(useInd)
            keepTimes{keepUID(j)} = tempTimes{useIndTemp(j)};
        end             
    end
    spikeT.UID = sort(keepUID);
    spikeT.times = keepTimes;
end

% only one type will be assigned to each cell
if ~isempty(type)
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
    keepUID = []; keepTimes = [];
    for i = 1:length(type)
        tempUID = []; tempTimes = [];
        tempUID = spikes.UID(contains(cell_metrics.putativeCellType, type(i)));
        tempTimes = spikes.times(contains(cell_metrics.putativeCellType, type(i)));
        [~,useIndTemp,useInd] = intersect(tempUID, spikeT.UID);
        keepUID = [keepUID spikeT.UID(useInd)];
        for j = 1:length(useInd)
            keepTimes{spikeT.UID(useInd(j))} = tempTimes{useIndTemp(j)};
        end            
    end
    spikeT.UID = sort(keepUID);
    spikeT.times = keepTimes;
end

% multiple states may be assigned to each cell's firing
if ~isempty(state)
    load(fullfile(basepath,[basename,'.SleepState.states.mat']));
    tempUID = spikeT.UID; tempTimes = spikeT.times; keepCt = 1;
    spikeT.UID = []; spikeT.times = [];
    for i = 1:length(tempUID)
        keepUn = false;
        for j = 1:length(state)
            if isfield(SleepState.ints, state(j))
                intervals = eval(strcat('SleepState.ints.',state(j)));
                tempSpk{i,j} = Restrict(tempTimes{tempUID(i)},intervals);
                if ~isempty(tempSpk{i,j})
                    keepUn = true;
                end
            else
                disp(fieldnames(SleepState.ints));
                error('Incorrect sleep state entered. Please choose from the above instead');
            end
        end
        if keepUn
            spikeT.UID(keepCt) = tempUID(i);
            spikeT.times{keepCt} = cat(1,tempSpk{i,:});
            spikeT.times{keepCt} = sort(spikeT.times{keepCt});
            keepCt = keepCt+1;
        end
    end        
end

% Might need to condense spikeT
if length(spikeT.times)>length(spikeT.UID)
    for i = 1:length(spikeT.UID)
        keepT{i} = spikeT.times{spikeT.UID(i)};
    end
    spikeT.times = []; spikeT.times = keepT;
end

end