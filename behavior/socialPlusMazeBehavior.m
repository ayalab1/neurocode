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
addParameter(p,'lapStart',1,@isnumeric); % percent of linear track to sep laps
addParameter(p,'speedTh',2,@isnumeric); % speed cm/sec threshold for stop/run times
addParameter(p,'savemat',true,@islogical); % save into animal.behavior.mat & linearTrackTrials.mat
addParameter(p,'show_fig',true,@islogical); % do you want a figure?
addParameter(p,'norm_zero_to_one',false,@islogical); % normalize linear coords 0-1
addParameter(p,'maze_sizes',[100 100],@isnumeric); % width of mazes in cm 
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
split_linearize = p.Results.split_linearize;
remove_extra_fields = p.Results.remove_extra_fields;
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

%% Conver to cm 
    if ~isempty(maze_sizes)
        pos_range = max(behavior.position.linearized) - min(behavior.position.linearized);
        convert_pix_to_cm_ratio = (pos_range / maze_sizes(1)); % using first maze size
        behavior.position.linearized = (behavior.position.linearized/ convert_pix_to_cm_ratio)';
        % convert xy to cm as well
        behavior.position.x = behavior.position.x /...
            convert_pix_to_cm_ratio;
        behavior.position.y = behavior.position.y /...
            convert_pix_to_cm_ratio;
        behavior.position.units = 'cm';
    end

    % normalize?
    if norm_zero_to_one
        behavior.position.linearized = ZeroToOne(behavior.position.linearized);
    end
   
%% Get laps
% needs to be done for each of the 4 arms independently 
behavior.posIndex = behavior.position.index;
behavior.posIndex(:,2) =NaN;
for a = 1:numel(unique(behavior.position.index))
    
    laps=FindLapsNSMAadapted(behavior.timestamps(behavior.position.index==a),...
        behavior.position.linearized(behavior.position.index==a),lapStart);

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
    
   dir1=find(ismember(behavior.timestamps',Restrict(behavior.timestamps(behavior.position.index==a),[inbound_start inbound_stop]))); 
   dir2=find(ismember(behavior.timestamps',Restrict(behavior.timestamps(behavior.position.index==a),[outbound_start outbound_stop]))); 
   behavior.posIndex(dir1,2)=1;
   behavior.posIndex(dir2,2)=2;
   
   % save lap information - NEEDS TO BE DONE DIFFERENTLY FOR THE 4 ARMS
    behavior.trials = [];
    behavior.trials =[[inbound_start;outbound_start],[inbound_stop;outbound_stop]];
    behavior.trials(1:length(inbound_start),3) = 1;
    behavior.trials(length(inbound_start)+1:length(inbound_start)+length(outbound_start),3) = 2;
    behavior.trials = sortrows(behavior.trials);
    behavior.trialID(1:length(behavior.trials),1) = a;
    behavior.trialID(:,2) = behavior.trials(:,3);
    behavior.trials = behavior.trials(:,1:2);
    behavior.trialIDname = {'centerToChamber';'chamberToCenter'}; % verify that this is correct 
    
    
end



% save speed threshold which might be different for each session
behavior.speedTh = speedTh;

%% Get periods of running
ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:)); 
t = behavior.timestamps(ok);
speed = LinearVelocity([behavior.timestamps(ok)', behavior.position.x(ok)', behavior.position.y(ok)'],5);
interpolated = Interpolate(speed,behavior.timestamps,'trim','off');
behavior.speed = interpolated(:,2)';

if isempty(maze_sizes)
    warning('data is not in cm, it is highly recommended to convert to cm')
    if isempty(speedTh)
        speedTh = 100;%prctile(behavior.speed(~isoutlier(behavior.speed)),10);
    end
end

run = t(FindInterval(interpolated(ok)>speedTh));
run = ConsolidateIntervals(run,'epsilon',0.01);
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
% remove outliers (data in between sessions gives outlier speeds)
[~,isOutlier] = RemoveOutliers(peak);
% remove run epochs that don't reach the speed threshold
run(peak<0.1 | isOutlier,:) = [];

% remove trial boundaries. The purpose of this is that the "turning" motion ending one trial would get separated from the
% running epoch on following trial and it can then be removed with the duration threshold
run = SubtractIntervals(run,bsxfun(@plus,sort(behavior.trials(:)),[-1 1]*0.001)); 
runDur = run(:,2)-run(:,1);
run(runDur<0.8 | runDur>15,:) = []; % remove run epochs that's too short or too long
behavior.run = run;

%% Separate positions for each direction of running
% this is the input that subsequent functions will use (e.g. findPlaceFieldsAvg1D)

% probably we don't need this
trials{1}.timestamps = [inbound_start inbound_stop];
%trials{1}.timestamps = trials{1}.timestamps(trials{1}.timestamps(:,2)-trials{1}.timestamps(:,1)<100,:); % excluding too long trials (need an input param)
trials{2}.timestamps = [outbound_start outbound_stop];
%trials{2}.timestamps = trials{2}.timestamps(trials{2}.timestamps(:,2)-trials{2}.timestamps(:,1)<100,:);

for i = 1:2
    behavior.positionTrials{i}=Restrict([behavior.timestamps' behavior.position.linearized'],trials{i}.timestamps);
    % probably we don't need this
    behavior.positionTrialsRun{i}=Restrict(behavior.positionTrials{i},run);
    trials{i}.timestampsRun = SubtractIntervals(trials{i}.timestamps, SubtractIntervals([0 Inf],run));
end

%% Manipulations


%% Plots to check results
% needs improvement
if show_fig
    figure;
    plot(behavior.timestamps,behavior.position.linearized-min(behavior.position.linearized),'k','LineWidth',2);hold on;
    plot(behavior.timestamps,behavior.speed,'r','LineWidth',2);hold on;
    PlotIntervals(trials{1}.timestamps,'color','b','alpha',.5);hold on;
    PlotIntervals(trials{2}.timestamps,'color','g','alpha',.5);hold on;
    ylim([-1,max(behavior.position.linearized)-min(behavior.position.linearized)])
    saveas(gcf,[basepath,filesep,[basename,'.linearTrackBehavior.fig']]);
end

if remove_extra_fields
    behavior = rmfield(behavior,{'positionTrials','run','positionTrialsRun'});
end

%% Generate output variables
if savemat
    if just_save_animal_behavior
        save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    else
        save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
        save([basepath,filesep,[basename,'.trials.mat']],'trials');
    end
end

end

