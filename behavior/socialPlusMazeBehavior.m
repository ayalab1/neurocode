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
%addParameter(p,'manipulation',[],@isstring); % add manipulation times to output
%addParameter(p,'lapStart',15,@isnumeric); % percent of linear track to sep laps
addParameter(p,'speedTh',5,@isnumeric); % speed cm/sec threshold for stop/run times
addParameter(p,'savemat',true,@islogical); % save into animal.behavior.mat & linearTrackTrials.mat
addParameter(p,'show_fig',true,@islogical); % do you want a figure?
addParameter(p,'norm_zero_to_one',false,@islogical); % normalize linear coords 0-1
addParameter(p,'maze_sizes',70,@isnumeric); % width of mazes in cm. 70 = calculated based on pythogorus theorem a2 = b2 + c2
addParameter(p,'just_save_animal_behavior',false,@islogical); % true will only save animal behav file

parse(p,varargin{:});
basepath = p.Results.basepath;
behavior = p.Results.behavior;
%manipulation = p.Results.manipulation;
%lapStart = p.Results.lapStart;
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

%% clean jumping tracker points/outliers (if any) from DeepLabCut
good_idx = manual_trackerjumps(behavior.timestamps,...
    behavior.position.x,...
    behavior.position.y,...
    startTime,...
    stopTime,...
    basepath,'darkmode',false);

behavior.position.x(~good_idx) = NaN;
behavior.position.y(~good_idx) = NaN;

%% Convert to cm
if ~isempty(behavior.maze_sizes)
    pos_range = max(behavior.position.x) - min(behavior.position.x);
    convert_pix_to_cm_ratio = (pos_range / behavior.maze_sizes(1)); % using first maze size
    
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
% 
% % normalize?
% if norm_zero_to_one
%     behavior.position.linearized = ZeroToOne(behavior.position.linearized);
% end

%% Find directions (centre to chamber = 1) and (chamber to center = -1) 
% find the difference between each position sucha that +ve value represents
% one direction and negative vale represents the other. 

linearized = behavior.position.linearized' ;% can apply median filter on linearized
%linearized = medfilt1(linearized)

linearizedDiff = diff(linearized);
linearizedAscend = linearizedDiff >= 0; 
linearizedAscendi = find(linearizedAscend) + 1; %add 1 to make a vector of same length
linearizedDescend = linearizedDiff < 0;
linearizedDescendi = find(linearizedDescend) + 1; %add 1 to make a vector of same length

dummy = nan(length(linearized),1); %create a dummy vector and add 1 or -1 for ascending and descendng tracker points.
dummy(linearizedAscendi) = 1;
dummy(linearizedDescendi) = -1;

% save direction and verify if makes sense
behavior.position.linearized_direction = dummy;
plot(behavior.timestamps(behavior.position.linearized_direction ==1),behavior.position.linearized(behavior.position.linearized_direction ==1),'.k');
hold on;
plot(behavior.timestamps(behavior.position.linearized_direction ==-1),behavior.position.linearized(behavior.position.linearized_direction ==-1),'.r');
hold off; 

%outbounds = behavior.position.linearized(dummy==1)
outboundi = dummy==1;
%inbounds= behavior.position.linearized(dummy==0)
inboundi = dummy==-1;

%verify if the ascending and descending tracking makes sense
plot(behavior.timestamps(outboundi),behavior.position.linearized(outboundi),'.k');
hold on;
plot(behavior.timestamps(inboundi),behavior.position.linearized(inboundi),'.r');
hold off; 

%% find timestamps intervals for outbounds(centre to chamber) and inbounds(chamber to center) 

intervalsOutbound = behavior.timestamps(FindInterval(outboundi));
%PlotIntervals(intervalsOutbound,'color','b','alpha',.5);hold on;
%sum(diff(intervalsOutbound,[],2))

%indexInterval = FindInterval(inboundi)
intervalsInbound = behavior.timestamps(FindInterval(inboundi));
%PlotIntervals(intervalsInbound,'color','b','alpha',.5);hold off;

%Create the variables to save in the general behavior file
behavior.direction = vertcat(intervalsOutbound,intervalsInbound);
behavior.directionLabel= [repmat(1,length(intervalsOutbound),1); repmat(-1,length(intervalsInbound),1)];

%%
% save speed threshold which might be different for each session
figure; plot(behavior.timestamps,behavior.speed, '.k');
PlotHVLines(behavior.speedTh,'h') % can increase or decrease speedTh depending on how well animal runs 

behavior.speedTh = speedTh;

%% Generate output variables
if savemat
    if just_save_animal_behavior
        save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    else
        save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    end
end

end
