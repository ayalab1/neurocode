function [firingMapsLaps] = firingMapLaps(positions,spikes,trials,varargin)

% USAGE
% [firingMaps] = firingMapLap(positions,spikes,trials,varargin)
% Calculates firing map for each lap for each neuron
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%   trials      - time intervals [trialstart trailstop] of each lap
%
%      or
%   behavior  - buzcode format behavior struct - NOT YET IMPLEMENTED
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold to compute firing rate
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: true) that saves firingMaps file
%     'CellInspector'  	- logical (default: false) that creates an otuput
%                   	compatible with CellInspector

%
%
% OUTPUT
%
%   firingMapLap - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%
% Can, Antonio FR, 3/2022

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'basepath',pwd,@ischar);

parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;
basepath = p.Results.basepath;

% number of conditions
conditions = numel(trials);
%%% TODO: conditions label

%% Calculate

% trials.timestampsRun contains running period already velocity filtered

%{
% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);
    
    % Compute timestamps where speed is under threshold
    positions{iCond}(v<speedThresh,:) = [];
end
%}

%%
% get firign rate maps
for unit = 1:length(spikes.times)
    for c = 1:conditions
        for t = 1:size(trials{1,1}.timestamps,1)
            
            trial_run_time = Restrict(trials{1,c}.timestampsRun,trials{1,c}.timestamps(t,:)); % find run time intervals of each trial
            position = Restrict(positions,trial_run_time); % find run position for each trial
            map{unit}{c}{t} = Map(position,spikes.times{unit},'smooth',smooth,'minTime',minTime,...
            'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
        end
    end
end
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMapsLaps.UID = spikes.UID;
try
    firingMapsLaps.sessionName = spikes.sessionName; % does not always exist in spikes
catch
    firingMapsLaps.sessionName = basenameFromBasepath(basepath);
end
try
    firingMapsLaps.region = spikes.region;
catch
    %warning('spikes.region is missing')
end

firingMapsLaps.params.smooth = smooth;
firingMapsLaps.params.minTime = minTime;
firingMapsLaps.params.nBins = nBins;
firingMapsLaps.params.x = map{1,1}{1,1}{1,1}.x;
firingMapsLaps.params.y = map{1,1}{1,1}{1,1}.y;
firingMapsLaps.params.maxGap = maxGap;
firingMapsLaps.params.mode = mode;
firingMapsLaps.params.maxDistance = maxDistance;

for unit = 1:length(spikes.times)
    for c = 1:conditions
        for t = 1:size(trials{1,1}.timestamps,1)
            firingMapsLaps.rateMaps{unit,1}{t,c} = map{unit}{c}{t}.z;
            firingMapsLaps.countMaps{unit,1}{t,c} = map{unit}{c}{t}.count;
            firingMapsLaps.occupancy{unit,1}{t,c} = map{unit}{c}{t}.time;
        end
    end
end

if saveMat
    save(fullfile(basepath,...
        [basenameFromBasepath(firingMapsLaps.sessionName),...
        '.firingMapsLaps.cellinfo.mat']),'firingMapsLaps');
end

end
