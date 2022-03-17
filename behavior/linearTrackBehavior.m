function [behavior,trials,posTrials] = linearTrackBehavior(varargin)
%        [behavior] = linearTrackBehavior(varargin)
% Gets raw tracking data and generates behavior structure
% based on the standards described in cellexplorer
% https://cellexplorer.org/datastructure/data-structure-and-format/#behavior

%  INPUTS  (Name-value paired inputs):
%
%  behavior = AYA lab standard behavior structure from general_behavior_file.m
%  tracking
%  manipulation = true if some type of manipulation (e.g. opto stimulation) was conducted
%  lapStart = start/end of lap, in % of track length
%  speedTh = speed threshold (times below it will be excluded)
%  savemat = save .mat varaibles to basepath
%
%  OUTPUTS
%
%  behavior = AYA lab standard behavior structure

% Can, Ryan, Antonio; 02/22
%
% TODO:
%       -find best way to store output trials. Ideally, should be easy to
%           identify, maybe table so human readable.
%       -Improve output figure
%       -add option to get laps from tracking or from sensors
%       -detect laps seperately per epoch in case linear track lengths are
%       different

%% parse inputs

p=inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'behavior',[],@isnumeric);
addParameter(p,'manipulation',true,@islogical); % add manipulation times to output
addParameter(p,'lapStart',20,@isnumeric); % percent of linear track to sep laps
addParameter(p,'speedTh',4,@isnumeric); % speed cm/sec threshold for stop/run times 
addParameter(p,'savemat',true,@islogical); % save into animal.behavior.mat & linearTrackTrials.mat
addParameter(p,'show_fig',true,@islogical); % do you want a figure?
addParameter(p,'norm_zero_to_one',false,@islogical); % normalize linear coords 0-1
addParameter(p,'maze_sizes',[],@isnumeric); % width of mazes in cm (must correspond with linear epochs)

parse(p,varargin{:});
basepath = p.Results.basepath;
behavior = p.Results.behavior;
manipulation = p.Results.manipulation;
lapStart = p.Results.lapStart;
speedTh = p.Results.speedTh;
savemat = p.Results.savemat;
show_fig = p.Results.show_fig;
norm_zero_to_one = p.Results.norm_zero_to_one;
maze_sizes = p.Results.maze_sizes;

basename = basenameFromBasepath(basepath);

%% Initialize/load behavior structure
% the basic animal.behavior.mat structure should have been generated first
% with general_behavior_file.m
if isempty(behavior)
    if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file')
        disp('detected animal.behavior.mat')
        load([basepath,filesep,[basename,'.animal.behavior.mat']]);
    else
        error('run general_behavior_file first')
    end
end

%% Pull in basename.session to epoch data
load([basepath,filesep,[basename,'.session.mat']]);
if ~isfield(session.epochs{1},'environment')
    warning('environment not labeled')
    warning('label environment and save before moving on')
    session = gui_session(session);
end
startTime = [];
stopTime = [];
for ep = session.epochs
    if contains(ep{1}.environment,'linear')
        startTime = [startTime;ep{1}.startTime];
        stopTime = [stopTime;ep{1}.stopTime];
    end
end
linear_epochs = [startTime,stopTime];

%% Linearize postions
% add method to select best LED/ marker

% make all coords outside linear track to nan
behavior.position.linearized = NaN(1,length(behavior.timestamps));

% make linear epoch by epoch to account for different maze positions
for ep = 1:size(linear_epochs,1)
    % find intervals in epoch
    [idx,~,~] = InIntervals(behavior.timestamps,linear_epochs(ep,:));
    % pull out xy
    xy = [behavior.position.x(idx)',behavior.position.y(idx)'];
    % make linear
    [~,lin,~] = pca(xy);
    linpos = lin(:,1);
    % min to zero
    linpos = linpos - min(linpos);
    if ~isempty(maze_sizes)
        pos_range = max(linpos) - min(linpos);
        convert_pix_to_cm_ratio = (pos_range / maze_sizes(ep));
        linpos = linpos / convert_pix_to_cm_ratio;
        % convert xy to cm as well
        behavior.position.x(idx) = behavior.position.x(idx) / convert_pix_to_cm_ratio;
        behavior.position.y(idx) = behavior.position.y(idx) / convert_pix_to_cm_ratio;
    end
    % normalize?
    if norm_zero_to_one
        linpos = ZeroToOne(linpos);
    end
    % add linear with interval epoch
    behavior.position.linearized(idx) = linpos';
end

%% Get laps
% add option to get laps from tracking or from sensors
laps=FindLapsNSMAadapted(behavior.timestamps,behavior.position.linearized,lapStart);

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

% save lap information 
behavior.trials = [];
behavior.trials =[[inbound_start;outbound_start],[inbound_stop;outbound_stop]];
behavior.trials(1:length(inbound_start),3) = 1;
behavior.trials(length(inbound_start)+1:length(inbound_start)+length(outbound_start),3) = 2;
behavior.trials = sortrows(behavior.trials);
behavior.trialID = behavior.trials(:,3);
behavior.trials = behavior.trials(:,1:2);
behavior.trialIDname = {'leftToRight';'rightToLeft'}; % verify that this is correct

%% Get periods of runnig
% this method works better than LinearVelocity.m
[~,~,~,vx,vy,~,~] = KalmanVel(behavior.position.linearized,behavior.position.linearized*0,behavior.timestamps,2);
v = sqrt(vx.^2+vy.^2); % TODO: v needs to be transformed to cm/s
behavior.speed =v'; % overriding the original speed variable (obtained with LinearVelocity.m)

if isempty(maze_sizes)
    warning('data is not standardized in cm')
    if isempty(speedTh)
        speedTh = prctile(v(~isoutlier(v)),10);
    end
end

% probably we don't need this
[quiet,quiescence] = QuietPeriods([behavior.timestamps' v],speedTh,0.5);
trials{1}.timestamps = [outbound_start outbound_stop];
trials{1}.timestamps = trials{1}.timestamps(trials{1}.timestamps(:,2)-trials{1}.timestamps(:,1)<100,:); % excluding too long trials (need an input param)
trials{2}.timestamps = [inbound_start inbound_stop];
trials{2}.timestamps = trials{2}.timestamps(trials{2}.timestamps(:,2)-trials{2}.timestamps(:,1)<100,:);
for i = 1:2
    trials{i}.timestampsRun = SubtractIntervals(trials{i}.timestamps,quiet);
end

%% Separate positions for each direction of running
% this is the input that subsequent functions will use (e.g. findPlaceFieldsAvg1D)
for i = 1:2
    behavior.positionTrials{i}=Restrict([behavior.timestamps' behavior.position.linearized'],trials{i}.timestamps);
    % probably we don't need this
    behavior.positionTrialsRun{i}=Restrict([behavior.timestamps' behavior.position.linearized'],trials{i}.timestampsRun);
end

%% Manipulations
if manipulation && exist([basepath,filesep,[basename,'.pulses.events.mat']],'file')
    load([basepath,filesep,[basename,'.pulses.events.mat']])
    behavior.manipulation = manipulation;
    behavior.stimON = pulses.intsPeriods;    
    % we need something more general because manipualtions can be stored in
    % digitalin or come in different ways
    
    % determine in which trials manipulations have been applied. For now I'm
    % only classifying each trial type (running directions in this case)
    % in stimON/OFF. We probably also want to do something like getting
    % the intersection of intervals for each trials with stimON
    stimTrials=[];
    for i = 1:numel(trials)
        t = InIntervals(behavior.stimON,trials{i}.timestampsRun);
        if sum(t) > 2
            stimTrials(i,1) = 1;
            behavior.trialIDname{i,2} = 'stimON';
            %trials{i}.manipulation = 'ON';
        else
            stimTrials(i,1) = 0;
            behavior.trialIDname{i,2} = 'stimOFF';
            %trials{i}.manipulation = 'OFF';
        end
        clear t;
    end
    for i = 1:numel(behavior.trialID)
        if behavior.trialID(i) == 1
           behavior.trialID(i,2) = stimTrials(1);
        elseif behavior.trialID(i) == 2
           behavior.trialID(i,2) = stimTrials(2);            
        end
    end
    
end

%% Plots to check results
% needs improvement
if show_fig
    figure;
    plot(behavior.timestamps,behavior.position.linearized,'k','LineWidth',2);hold on;
    plot(behavior.timestamps,v,'r','LineWidth',2);hold on;
    PlotIntervals(behavior.trials(behavior.trialID(:,1)==1,:),'color','b','alpha',.5);hold on;
    PlotIntervals(behavior.trials(behavior.trialID(:,1)==2,:),'color','g','alpha',.5);hold on;
    if ~isempty(manipulation) && exist([basepath,filesep,[basename,'.pulses.events.mat']],'file')
        PlotIntervals(pulses.intsPeriods,'color','k','alpha',.5);hold on;
    end
    ylim([min(behavior.position.linearized),max(behavior.position.linearized)])
    saveas(gcf,[basepath,filesep,[basename,'.linearTrackBehavior.fig']]);
end

%% Generate output variables
if savemat
    save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
end

end

