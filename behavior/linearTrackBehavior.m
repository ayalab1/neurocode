function [behavior,trials,posTrials] = linearTrackBehavior(varargin)
%        [behavior] = linearTrackBehavior(varargin)
% Gets raw tracking data and generates behavior structure
% based on the standards described in cellexplorer
% https://cellexplorer.org/datastructure/data-structure-and-format/#behavior

%  INPUTS  (Name-value paired inputs):
%
%  tracking
%  manipulation = true if some type of manipulation (e.g. opto stimulation) was conducted
%  lapStart = start/end of lap, in % of track length
%  speedTh = speed threshold (times below it will be excluded)
%  savemat = save .mat varaibles to basepath
%
%  OUTPUTS
%
%  behavior = AYA lab standard behavior structure
%  trials and postTrials are just temporary, should not be needed once this
%  function is finished

% Can, Ryan, Antonio; 02/22
%
% TODO: 
%       -data should be restricted to just linear track epochs
%       -xy should be converted from pixels to cm (speed thres can them be adjusted)
%       -find best way to store output trials. Ideally, should be easy to
%           identify, maybe table so human readable.
%       -Improve output figure
%       -add option to get laps from tracking or from sensors
%
%% parse inputs

p=inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'tracking',[],@isnumeric);
addParameter(p,'manipulation',true,@islogical);
addParameter(p,'lapStart',20,@isnumeric);
addParameter(p,'speedTh',0.1,@isnumeric);
addParameter(p,'savemat',true,@islogical);
addParameter(p,'show_fig',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
tracking = p.Results.tracking;
manipulation = p.Results.manipulation;
lapStart = p.Results.lapStart;
speedTh = p.Results.speedTh;
savemat = p.Results.savemat;
show_fig = p.Results.show_fig;

basename = basenameFromBasepath(basepath);

%% Initialize/load behavior structure
% not sure what do we prefer to do here. We could generate the behavior
% structure with Ryan's function and just load it here or perhaps better
% generate it directly here from Tracking.behavior.mat
if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file')
    disp('detected animal.behavior.mat')
    load([basepath,filesep,[basename,'.animal.behavior.mat']]);
else
    error('run general_behavior_file first')
end

%% Linearize postions
% add method to select best LED/ marker
[~,lin,~] = pca([behavior.position.x',behavior.position.y']);
linpos = lin(:,1);
linpos = ZeroToOne(linpos); % TODO xy needs to be transformed to cm
behavior.position.linearized = linpos';

%% Get laps
% add option to get laps from tracking or from sensors
laps=FindLapsNSMAadapted(behavior.time,linpos,lapStart);

outbound_start=[];outbound_stop=[];inbound_start=[];inbound_stop=[];
for i = 1:length([laps.start_ts])-1
    if laps(i).direction == 1
        outbound_start=cat(1,outbound_start,laps(i).start_ts);
        outbound_stop=cat(1,outbound_stop,laps(i+1).start_ts);
    elseif laps(i).direction == -1
        inbound_start=cat(1,inbound_start,laps(i).start_ts);
        inbound_stop=cat(1,inbound_stop,laps(i+1).start_ts);
    end
end

% I am creating for now a temporary trials struct, but this needs to be
% included inside the behavior struct instead
trials{1}.timestamps = [outbound_start outbound_stop];
trials{1}.timestamps = trials{1}.timestamps(trials{1}.timestamps(:,2)-trials{1}.timestamps(:,1)<100,:); % excluding too long trials (need an input param)
trials{2}.timestamps = [inbound_start inbound_stop];
trials{2}.timestamps = trials{2}.timestamps(trials{2}.timestamps(:,2)-trials{2}.timestamps(:,1)<100,:);

%% Get periods of runnig
% this method works better than LinearVelocity.m
[~,~,~,vx,vy,~,~] = KalmanVel(linpos,linpos*0,behavior.time,2);
v = sqrt(vx.^2+vy.^2); % TODO: v needs to be transformed to cm/s

[quiet,quiescence] = QuietPeriods([behavior.time' v],speedTh,0.5);

for i = 1:2
    trials{i}.timestampsRun = SubtractIntervals(trials{i}.timestamps,quiet);
end

%% Separate positions for each direction of running
% this is the input that subsequent functions will use (e.g. findPlaceFieldsAvg1D)
% so, for now, I am converting it to the old posTrials. The best solution
% would be to modify findPlaceFieldsAvg1D and other functions to take
% instead the behavior structure
for i = 1:2
    trials{i}.positions=Restrict([behavior.time' linpos],trials{i}.timestamps);
    trials{i}.positionsRun=Restrict([behavior.time' linpos],trials{i}.timestampsRun);
    posTrials{i} = trials{i}.positions;
end

%% Manipulations
if manipulation && exist([basepath,filesep,[basename,'.pulses.events.mat']],'file')
    load([basepath,filesep,[basename,'.pulses.events.mat']])
    stimON = pulses.intsPeriods;
    % we need something more general because manipualtions can be stored in
    % digitalin or come in different ways
    
    % determine in which trials manipulations have been applied. For now I'm
    % only classifying each trial type (running directions in this case)
    % in stimON/OFF. We probably also want to do something like getting
    % the intersection of intervals for each trials with stimON
    for i = 1:numel(trials)
        t = InIntervals(stimON,trials{i}.timestampsRun);
        if sum(t) > 1
            trials{i}.manipulation = 'ON';
        else
            trials{i}.manipulation = 'OFF';
        end
        clear t;
    end
end

%% Plots to check results
% need improvement
if show_fig
    figure;
    plot(behavior.time,behavior.position.linearized,'k','LineWidth',2);hold on;
    plot(behavior.time,v,'r','LineWidth',2);hold on;
    PlotIntervals(trials{1}.timestampsRun,'color','b','alpha',.5);hold on;
    PlotIntervals(trials{2}.timestampsRun,'color','g','alpha',.5);hold on;
    if manipulation && exist([basepath,filesep,[basename,'.pulses.events.mat']],'file')
        PlotIntervals(pulses.intsPeriods,'color','m','alpha',.5);hold on;
    end
    ylim([0 1]);

    saveas(gcf,[basepath,filesep,[basename,'.linearTrackBehavior.fig']]);
end
%% Generate output variables
if savemat
    save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    save([basepath,filesep,[basename,'.linearTrackTrials.mat']],'trials','posTrials');
end

end

