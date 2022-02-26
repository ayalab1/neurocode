function [behavior,trials,posTrials] = linearTrackBehavior(varargin)
%        [behavior] = linearTrackBehavior(varargin)
% Gets raw tracking data and generates behavior structure
% based on the standards described in cellexplorer
% https://cellexplorer.org/datastructure/data-structure-and-format/#behavior

%  INPUTS  (Name-value paired inputs):

%  tracking
%  lapStart = start/end of lap, in % of track length
%  speedTh = speed threshold (times below it will be excluded)

% Can, Ryan, Antonio; 02/22

%% parse inputs

p=inputParser;
addParameter(p,'basepath',pwd);
addParameter(p,'tracking',[],@isnumeric);
addParameter(p,'lapStart',20,@isnumeric);
addParameter(p,'speedTh',0.1,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
tracking = p.Results.tracking;
lapStart = p.Results.lapStart;
speedTh = p.Results.speedTh;

cd(basepath)
basename = basenameFromBasepath(pwd);

if isempty(tracking)
    load([basename '.Tracking.behavior.mat']);
    % why there is Tracking.bahavior and tracking.behavior ??
end

%% Initialize/load behavior structure
% not sure what do we prefer to do here. We could generate the behavior
% structure with Ryan's function and just load it here or perhaps better
% generate it directly here from Tracking.behavior.mat


if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file')
    disp('detected animal.behavior.mat')
    load([basepath,filesep,[basename,'.animal.behavior.mat']])
end

%% Linearize postions 
% add a method to select best LED/ marker 
[~,lin,~] = pca([behavior.position.x',behavior.position.y']);
linpos = lin(:,1); linpos = ZeroToOne(linpos);
behavior.linearized = linpos;

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
% this method is better than LinearVelocity.m
[~,~,~,vx,vy,~,~] = KalmanVel(linpos,linpos*0,behavior.time,2);
v = sqrt(vx.^2+vy.^2);  

[quiet,quiescence] = QuietPeriods([behavior.time' v],0.1,0.5);

for i = 1:2
    trials{i}.timestampsRun = SubtractIntervals(trials{i}.timestamps,quiet);
end

%% Separate positions for each direction of running
% this is the input that subsequent functions will use (e.g. findPlaceFieldsAvg1D)
% so, for now, I am converting it to the old posTrials 
for i = 1:2
    trials{i}.positions=Restrict([behavior.time' linpos],trials{i}.timestamps);
    trials{i}.positionsRun=Restrict([behavior.time' linpos],trials{i}.timestampsRun); 
    posTrials{i} = trials{i}.positions;
end
   
%% Plots to check results 
% need improvement
figure;plot(behavior.time,behavior.linearized,'k','LineWidth',2);hold on;
       plot(behavior.time,v,'r','LineWidth',2);hold on;
PlotIntervals(trials{1}.timestampsRun,'color','b','alpha',.5);hold on;
PlotIntervals(trials{2}.timestampsRun,'color','g','alpha',.5);hold on;

% testing stuff...
if exist([basepath,filesep,[basename,'.pulses.events.mat']],'file')
    load([basepath,filesep,[basename,'.pulses.events.mat']])
end
PlotIntervals(pulses.intsPeriods,'color','m','alpha',.5);hold on;

%% Generate output variables 


end

