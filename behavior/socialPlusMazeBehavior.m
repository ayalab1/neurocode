function [behavior] = socialPlusMazeBehavior(behavior,varargin)
%        [behavior] = socialPlusMazeBehavior(behavior,varargin)
% Takes basic behavior structure generated with general_behavior_file.m and
% completes it for the plus maze social behavior task. Assumes that
% behavior.mat already contains linearrized positions separated for the
% four maze arms

%  INPUTS  (Name-value paired inputs):
%
%  behavior = AYA lab standard behavior structure from general_behavior_file.m
%  manipulation = true if some type of manipulation (e.g. opto stimulation) was conducted
%  lapStart = start/end of lap, in % of track length
%  speedTh = speed threshold (times below it will be excluded)
%  savemat = save .mat varaibles to basepath
%
%  OUTPUTS
%
%  behavior = AYA lab standard behavior structure

%  Antonio, Praveen 07/22
%
% TODO:


%% parse inputs

p=inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'behavior',[],@isnumeric);
addParameter(p,'manipulation',[],@isstring); % add manipulation times to output
addParameter(p,'lapStart',15,@isnumeric); % percent of linear track to sep laps
addParameter(p,'speedTh',2,@isnumeric); % speed cm/sec threshold for stop/run times
addParameter(p,'savemat',true,@islogical); % save into animal.behavior.mat & linearTrackTrials.mat
addParameter(p,'show_fig',true,@islogical); % do you want a figure?
addParameter(p,'norm_zero_to_one',false,@islogical); % normalize linear coords 0-1
addParameter(p,'maze_sizes',70,@isnumeric); % width of mazes in cm. 70 = calculated based on pythogorus theorem a2 = b2 + c2
addParameter(p,'just_save_animal_behavior',false,@islogical); % true will only save animal behav file

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
%split_linearize = p.Results.split_linearize;
%remove_extra_fields = p.Results.remove_extra_fields;
just_save_animal_behavior = p.Results.just_save_animal_behavior;

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
    if contains(ep{1}.environment,'social_plusmaze')
        startTime = [startTime;ep{1}.startTime];
        stopTime = [stopTime;ep{1}.stopTime];
    end
end
task_epochs = [startTime,stopTime];

%% Convert to cm
if ~isempty(maze_sizes)
    pos_range = max(behavior.position.x) - min(behavior.position.x);
    convert_pix_to_cm_ratio = (pos_range / maze_sizes(1)); % using first maze size
    
    behavior.position.linearized = (behavior.position.linearized / convert_pix_to_cm_ratio)';
    
    % convert xy to cm as well
    behavior.position.x = behavior.position.x /...
        convert_pix_to_cm_ratio;
    behavior.position.y = behavior.position.y /...
        convert_pix_to_cm_ratio;
    
    behavior.position.units = 'cm';
    
    % convert all other points to cm
    fields_labels = fields(behavior.position);
    
    for label = fields_labels(contains(fields(behavior.position),'_x_'))'
        behavior.position.(label{1}) = behavior.position.(label{1}) /...
            convert_pix_to_cm_ratio;
    end
    for label = fields_labels(contains(fields(behavior.position),'_y_'))'
        behavior.position.(label{1}) = behavior.position.(label{1}) /...
            convert_pix_to_cm_ratio;
    end
end

% normalize?
if norm_zero_to_one
    behavior.position.linearized = ZeroToOne(behavior.position.linearized);
end

%% Get laps
states = unique(behavior.position.index(~isnan(behavior.position.index)));

%save states and tracks IDs/names in behavior file
behavior.states = {0,'trackA';1,'trackB'; 2,'trackC'; 3,'trackD'};
behavior.trialsIDname = {1,'chamberToCenter';2, 'centerToChamber'};


for a = states(1)' % trials for trackA
    
    laps = FindLapsNSMAadapted(behavior.timestamps(behavior.position.index==a),...
        behavior.position.linearized(behavior.position.index==a),lapStart);

    outbound_start=[]; outbound_stop=[]; inbound_start=[]; inbound_stop=[];
    for i = 1:length([laps.start_ts])-1
        if laps(i).direction == 1
            outbound_start=cat(1,outbound_start,laps(i).start_ts);
            outbound_stop=cat(1,outbound_stop,laps(i+1).start_ts);
        elseif laps(i).direction == -1
            inbound_start=cat(1,inbound_start,laps(i).start_ts);
            inbound_stop=cat(1,inbound_stop,laps(i+1).start_ts);
        end
    end
    
    behavior.trials_trackA = [];
    behavior.trials_trackA =[[inbound_start;outbound_start],[inbound_stop;outbound_stop]];
    behavior.trials_trackA(1:length(inbound_start),3) = 1;
    behavior.trials_trackA(length(inbound_start)+1:length(inbound_start)+length(outbound_start),3) = 2;
    behavior.trials_trackA = sortrows(behavior.trials_trackA,1);

    % save trials information in behavior file
    behavior.trials_trackA(:,3) = behavior.trials_trackA(:,3);
    behavior.trials_trackA(:,1:2) = behavior.trials_trackA(:,1:2);
end

for a = states(2)' % trials for trackB
    
    laps = FindLapsNSMAadapted(behavior.timestamps(behavior.position.index==a),...
        behavior.position.linearized(behavior.position.index==a),lapStart);

    outbound_start=[]; outbound_stop=[]; inbound_start=[]; inbound_stop=[];
    for i = 1:length([laps.start_ts])-1
        if laps(i).direction == 1
            outbound_start=cat(1,outbound_start,laps(i).start_ts);
            outbound_stop=cat(1,outbound_stop,laps(i+1).start_ts);
        elseif laps(i).direction == -1
            inbound_start=cat(1,inbound_start,laps(i).start_ts);
            inbound_stop=cat(1,inbound_stop,laps(i+1).start_ts);
        end
    end
    
    behavior.trials_trackB = [];
    behavior.trials_trackB =[[inbound_start;outbound_start],[inbound_stop;outbound_stop]];
    behavior.trials_trackB(1:length(inbound_start),3) = 1;
    behavior.trials_trackB(length(inbound_start)+1:length(inbound_start)+length(outbound_start),3) = 2;
    behavior.trials_trackB = sortrows(behavior.trials_trackB,1);

    % save trials information in behavior file
    behavior.trials_trackB(:,3) = behavior.trials_trackB(:,3);
    behavior.trials_trackB(:,1:2) = behavior.trials_trackB(:,1:2);
end

for a = states(3)' % trials for trackC
    
    laps = FindLapsNSMAadapted(behavior.timestamps(behavior.position.index==a),...
        behavior.position.linearized(behavior.position.index==a),lapStart);

    outbound_start=[]; outbound_stop=[]; inbound_start=[]; inbound_stop=[];
    for i = 1:length([laps.start_ts])-1
        if laps(i).direction == 1
            outbound_start=cat(1,outbound_start,laps(i).start_ts);
            outbound_stop=cat(1,outbound_stop,laps(i+1).start_ts);
        elseif laps(i).direction == -1
            inbound_start=cat(1,inbound_start,laps(i).start_ts);
            inbound_stop=cat(1,inbound_stop,laps(i+1).start_ts);
        end
    end
    
    behavior.trials_trackC = [];
    behavior.trials_trackC =[[inbound_start;outbound_start],[inbound_stop;outbound_stop]];
    behavior.trials_trackC(1:length(inbound_start),3) = 1;
    behavior.trials_trackC(length(inbound_start)+1:length(inbound_start)+length(outbound_start),3) = 2;
    behavior.trials_trackC = sortrows(behavior.trials_trackC,1);

    % save trials information in behavior file
    behavior.trials_trackC(:,3) = behavior.trials_trackC(:,3);
    behavior.trials_trackC(:,1:2) = behavior.trials_trackC(:,1:2);
end

for a = states(4)' % trials for trackD
    
    laps = FindLapsNSMAadapted(behavior.timestamps(behavior.position.index==a),...
        behavior.position.linearized(behavior.position.index==a),lapStart);

    outbound_start=[]; outbound_stop=[]; inbound_start=[]; inbound_stop=[];
    for i = 1:length([laps.start_ts])-1
        if laps(i).direction == 1
            outbound_start=cat(1,outbound_start,laps(i).start_ts);
            outbound_stop=cat(1,outbound_stop,laps(i+1).start_ts);
        elseif laps(i).direction == -1
            inbound_start=cat(1,inbound_start,laps(i).start_ts);
            inbound_stop=cat(1,inbound_stop,laps(i+1).start_ts);
        end
    end
    
    behavior.trials_trackD = [];
    behavior.trials_trackD =[[inbound_start;outbound_start],[inbound_stop;outbound_stop]];
    behavior.trials_trackD(1:length(inbound_start),3) = 1;
    behavior.trials_trackD(length(inbound_start)+1:length(inbound_start)+length(outbound_start),3) = 2;
    behavior.trials_trackD = sortrows(behavior.trials_trackD,1);

    % save trials information in behavior file
    behavior.trials_trackD(:,3) = behavior.trials_trackD(:,3);
    behavior.trials_trackD(:,1:2) = behavior.trials_trackD(:,1:2);
end

% save speed threshold which might be different for each session
behavior.speedTh = speedTh;

%% Get periods of running
ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:));
t = behavior.timestamps(ok);
speed = LinearVelocity([behavior.timestamps(ok)', behavior.position.x(ok)', behavior.position.y(ok)'],5);
interpolated = Interpolate(speed,behavior.timestamps,'trim','off');
behavior.speed = interpolated(:,2)';

% if isempty(maze_sizes)
%     warning('data is not in cm, it is highly recommended to convert to cm')
%     if isempty(speedTh)
%         speedTh = 100;%prctile(behavior.speed(~isoutlier(behavior.speed)),10);
%     end
% end

% run = t(FindInterval(interpolated(ok)>speedTh));
% run = ConsolidateIntervals(run,'epsilon',0.01);
% [in,w] = InIntervals(behavior.timestamps(:),run);
% peak = Accumulate(w(in),behavior.speed(in)','mode','max');
% % remove outliers (data in between sessions gives outlier speeds)
% [~,isOutlier] = RemoveOutliers(peak);
% % remove run epochs that don't reach the speed threshold
% run(peak<0.1 | isOutlier,:) = [];
% 
% % remove trial boundaries. The purpose of this is that the "turning" motion ending one trial would get separated from the
% % running epoch on following trial and it can then be removed with the duration threshold
% run = SubtractIntervals(run,bsxfun(@plus,sort(behavior.trials(:)),[-1 1]*0.001));
% runDur = run(:,2)-run(:,1);
% run(runDur<0.8 | runDur>15,:) = []; % remove run epochs that's too short or too long
% behavior.run = run;

%% Manipulations


%% Plots to check results
%NEED TO IMPLEMENT
if show_fig
    figure;
    plot(behavior.timestamps,behavior.position.linearized-min(behavior.position.linearized),'.k','LineWidth',2);hold on;
    plot(behavior.timestamps,behavior.speed,'r','LineWidth',2);hold on;
    PlotIntervals(behavior.trials(behavior.trialID == 1,:),'color','b','alpha',.5);hold on;
    PlotIntervals(behavior.trials(behavior.trialID == 2,:),'color','g','alpha',.5);hold on;
    ylim([-1,max(behavior.position.linearized)-min(behavior.position.linearized)])
    saveas(gcf,[basepath,filesep,[basename,'.linearTrackBehavior.fig']]);
end

% if remove_extra_fields
%     behavior = rmfield(behavior,{'positionTrials','run','positionTrialsRun'});
% end

%% Generate output variables
if savemat
    if just_save_animal_behavior
        save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    else
        save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
%         save([basepath,filesep,[basename,'.trials.mat']],'trials');
    end
end

end
