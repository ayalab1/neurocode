function [behavior] = linearTrackBehavior(varargin)

% Generate behavioral file for linear track sessions.
%
% [behavior] = linearTrackBehavior(varargin)
% Gets tracking data and adds to behavior structure
% based on the standards described in CellExplorer
% https://cellexplorer.org/datastructure/data-structure-and-format/#behavior

%  INPUTS  (Name-value paired inputs):
%
%  behavior = AYA lab standard behavior structure from general_behavior_file.m
%  tracking
%  manipulation = true if some type of manipulation (e.g. opto stimulation)
%  was conducted  - (Needs better explanation, the way it is written makes
%  confusing how to implement- HLR 2/24)
%  lapStart = start/end of lap, in % of track length
%  speedTh = speed threshold (times below it will be excluded)
%  savemat = save .mat varaibles to basepath
%
%  OUTPUTS
%
%  behavior = cellexplorer standard behavior structure
%
% % Copyright (C) 2022 Can Liu, Ryan Harvey, Antonio FR
%
% TODO:
%       -Current implementation does not entirely conform to cell explorer:
%           https://cellexplorer.org/datastructure/data-structure-and-format/#behavior
%       -Improve output figure
%       -add option to get laps from tracking or from sensors

%% parse inputs

p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder);
addParameter(p, 'behavior', [], @isnumeric);
addParameter(p, 'manipulation', [], @isstring); % add manipulation times to output
addParameter(p, 'lapStart', 20, @isnumeric); % percent of linear track to sep laps
addParameter(p, 'speedTh', 4, @isnumeric); % speed cm/sec threshold for stop/run times
addParameter(p, 'savemat', true, @islogical); % save into animal.behavior.mat & linearTrackTrials.mat
addParameter(p, 'show_fig', true, @islogical); % do you want a figure?
addParameter(p, 'norm_zero_to_one', false, @islogical); % normalize linear coords 0-1
addParameter(p, 'maze_sizes', [], @isnumeric); % width of mazes in cm (must correspond with linear epochs)
addParameter(p, 'split_linearize', false, @islogical); % make linear epoch by epoch
addParameter(p, 'remove_extra_fields', false, @islogical); % removes extra FMA syle fields 'positionTrials','run','positionTrialsRun'
addParameter(p, 'just_save_animal_behavior', false, @islogical); % true will only save animal behav file
addParameter(p, 'clean_tracker_jumps', false, @islogical);

parse(p, varargin{:});
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
clean_tracker_jumps = p.Results.clean_tracker_jumps;

basename = basenameFromBasepath(basepath);

%% Initialize/load behavior structure
% the basic animal.behavior.mat structure should have been generated first
% with general_behavior_file.m
if isempty(behavior)
    if exist([basepath, filesep, [basename, '.animal.behavior.mat']], 'file')
        disp('detected animal.behavior.mat')
        load([basepath, filesep, [basename, '.animal.behavior.mat']], 'behavior');
    else
        error('run general_behavior_file first')
    end
end

%% Pull in basename.session to epoch data
load([basepath, filesep, [basename, '.session.mat']], 'session');
if ~isfield(session.epochs{1}, 'environment')
    warning('environment not labeled')
    warning('label environment and save before moving on')
    session = gui_session(session);
end
startTime = [];
stopTime = [];
for ep = session.epochs
    if contains(lower(ep{1}.environment), 'linear')
        startTime = [startTime; ep{1}.startTime];
        stopTime = [stopTime; ep{1}.stopTime];
    end
end
linear_epochs = [startTime, stopTime];

%% Linearize postions
% make all coords outside linear track to nan
if isempty(behavior.position.linearized)
    behavior.position.linearized = NaN(1, length(behavior.position.x));
end

% make linear all linear track coords together across epochs
if ~split_linearize
    xy = [];
    idxs = [];
    for ep = 1:size(linear_epochs, 1)
        % find intervals in epoch
        [idx, ~, ~] = InIntervals(behavior.timestamps, linear_epochs(ep, :));
        % pull out xy
        xy = [xy; behavior.position.x(idx)', behavior.position.y(idx)'];
        % store index
        idxs = [idxs, idx];
    end
    idxs = sum(idxs, 2) > 0;

    % make linear
    [~, lin, ~] = pca(xy);
    linpos = lin(:, 1);
    % min to zero
    %linpos = linpos - min(linpos);
    if ~isempty(maze_sizes)
        pos_range = max(linpos) - min(linpos);
        convert_pix_to_cm_ratio = (pos_range / maze_sizes(1)); % using first maze size
        linpos = linpos / convert_pix_to_cm_ratio;
        % convert xy to cm as well
        behavior.position.x(idxs) = behavior.position.x(idxs) / ...
            convert_pix_to_cm_ratio;
        behavior.position.y(idxs) = behavior.position.y(idxs) / ...
            convert_pix_to_cm_ratio;
        behavior.position.units = 'cm';
    end
    % normalize?
    if norm_zero_to_one
        linpos = ZeroToOne(linpos);
    end
    % add linear with interval epoch
    behavior.position.linearized(idxs) = linpos';

else % if want to linearize tracking epoch by epoch

    % make linear epoch by epoch to account for different maze positions
    for ep = 1:size(linear_epochs, 1)
        % find intervals in epoch
        [idx, ~, ~] = InIntervals(behavior.timestamps, linear_epochs(ep, :));
        % pull out xy
        xy = [behavior.position.x(idx)', behavior.position.y(idx)'];
        % make linear
        [~, lin, ~] = pca(xy);
        linpos = lin(:, 1);
        % min to zero
        %linpos = linpos - min(linpos);
        if ~isempty(maze_sizes)
            pos_range = max(linpos) - min(linpos);
            convert_pix_to_cm_ratio = (pos_range / maze_sizes(ep));
            linpos = linpos / convert_pix_to_cm_ratio;
            % convert xy to cm as well
            behavior.position.x(idx) = behavior.position.x(idx) / ...
                convert_pix_to_cm_ratio;
            behavior.position.y(idx) = behavior.position.y(idx) / ...
                convert_pix_to_cm_ratio;
            behavior.position.units = 'cm';
        end
        % normalize?
        if norm_zero_to_one
            linpos = ZeroToOne(linpos);
        end
        % add linear with interval epoch
        behavior.position.linearized(idx) = linpos';
    end
end

% Clean:
if clean_tracker_jumps
    xy = ([behavior.position.x(:), behavior.position.y(:)]);
    smoothed = nansmooth(xy, [10, 0]);
    d = sqrt(sum((xy - smoothed).^2, 2));
    [~, ~, threshold] = isoutlier(d, 'quartiles', 'ThresholdFactor', 5);
    bad = isnan(xy(:, 1)) | (d > threshold);
    behavior.position.x(bad) = nan;
    behavior.position.y(bad) = nan;
    behavior.position.linearized(bad) = nan;
end

%% Get laps
% add option to get laps from tracking or from sensors
if ~split_linearize
    laps = FindLapsNSMAadapted(behavior.timestamps, behavior.position.linearized, lapStart);

    outbound_start = [];
    outbound_stop = [];
    inbound_start = [];
    inbound_stop = [];
    for i = 1:length([laps.start_ts]) - 1
        if laps(i).direction == 1
            outbound_start = cat(1, outbound_start, laps(i).start_ts);
            outbound_stop = cat(1, outbound_stop, laps(i+1).start_ts);
        elseif laps(i).direction == -1
            inbound_start = cat(1, inbound_start, laps(i).start_ts);
            inbound_stop = cat(1, inbound_stop, laps(i+1).start_ts);
        end
    end
    inbound_intervals = [inbound_start, inbound_stop];
    outbound_intervals = [outbound_start, outbound_stop];
else
    inbound_intervals = [];
    outbound_intervals = [];
    for ep = 1:size(linear_epochs, 1)
        % find intervals in epoch
        [idx, ~, ~] = InIntervals(behavior.timestamps, linear_epochs(ep, :));
        laps = FindLapsNSMAadapted(behavior.timestamps(idx), behavior.position.linearized(idx), lapStart);

        outbound_start = [];
        outbound_stop = [];
        inbound_start = [];
        inbound_stop = [];
        for i = 1:length([laps.start_ts]) - 1
            if laps(i).direction == 1
                outbound_start = cat(1, outbound_start, laps(i).start_ts);
                outbound_stop = cat(1, outbound_stop, laps(i+1).start_ts);
            elseif laps(i).direction == -1
                inbound_start = cat(1, inbound_start, laps(i).start_ts);
                inbound_stop = cat(1, inbound_stop, laps(i+1).start_ts);
            end
        end
        inbound_intervals = [inbound_intervals; [inbound_start, inbound_stop]];
        outbound_intervals = [outbound_intervals; [outbound_start, outbound_stop]];
    end
end

% save lap information
behavior.trials = [inbound_intervals; outbound_intervals];
behavior.trials(1:size(inbound_intervals, 1), 3) = 1;
behavior.trials(size(inbound_intervals, 1)+1:size(inbound_intervals, 1)+size(outbound_intervals, 1), 3) = 2;
behavior.trials = sortrows(behavior.trials);
behavior.trialID = behavior.trials(:, 3);
behavior.trials = behavior.trials(:, 1:2);
behavior.trialIDname = {'leftToRight'; 'rightToLeft'}; % verify that this is correct

% save speed threshold which might be different for each session

behavior.speedTh = speedTh;

%% Get periods of running
speed = LinearVelocity([behavior.timestamps', behavior.position.x', behavior.position.y'], 5);
behavior.speed = speed(:, 2)';

if isempty(maze_sizes)
    if isempty(speedTh)
        speedTh = 100; %prctile(behavior.speed(~isoutlier(behavior.speed)),10);
    end
end

run = behavior.timestamps(FindInterval(behavior.speed > speedTh));
run = ConsolidateIntervals(run, 'epsilon', 0.01);
[in, w] = InIntervals(behavior.timestamps(:), run);
peak = Accumulate(w(in), behavior.speed(in)', 'mode', 'max');
% remove outliers (data in between sessions gives outlier speeds)
[~, isOutlier] = RemoveOutliers(peak, 10); % set the "outlier" threshold pretty high to make sure high-speed running epochs don't get removed
% remove run epochs that don't reach the speed threshold
run(peak < 0.1 | isOutlier, :) = [];

% remove trial boundaries. The purpose of this is that the "turning" motion ending one trial would get separated from the
% running epoch on following trial and it can then be removed with the duration threshold
run = SubtractIntervals(run, bsxfun(@plus, sort(behavior.trials(:)), [-1, 1]*0.001));
runDur = run(:, 2) - run(:, 1);
run(runDur < 0.8 | runDur > 15, :) = []; % remove run epochs that's too short or too long
behavior.run = run;

%% Separate positions for each direction of running
% this is the input that subsequent functions will use (e.g. findPlaceFieldsAvg1D)

% probably we don't need this
trials{1}.timestamps = inbound_intervals;
%trials{1}.timestamps = trials{1}.timestamps(trials{1}.timestamps(:,2)-trials{1}.timestamps(:,1)<100,:); % excluding too long trials (need an input param)
trials{2}.timestamps = outbound_intervals;
%trials{2}.timestamps = trials{2}.timestamps(trials{2}.timestamps(:,2)-trials{2}.timestamps(:,1)<100,:);

for i = 1:2
    behavior.positionTrials{i} = Restrict([behavior.timestamps', behavior.position.linearized'], trials{i}.timestamps);
    % probably we don't need this
    behavior.positionTrialsRun{i} = Restrict(behavior.positionTrials{i}, run);
    trials{i}.timestampsRun = SubtractIntervals(trials{i}.timestamps, SubtractIntervals([0, Inf], run));
end

%% Manipulations
if ~isempty(manipulation) && exist([basepath, filesep, [basename, '.pulses.events.mat']], 'file')
    load([basepath, filesep, [basename, '.pulses.events.mat']], 'pulses')
    behavior.manipulation = manipulation;
    behavior.stimON = pulses.intsPeriods;
    % we need something more general because manipualtions can be stored in
    % digitalin or come in different ways

    % determine in which trials manipulations have been applied. For now I'm
    % only classifying each trial type (running directions in this case)
    % in stimON/OFF. We probably also want to do something like getting
    % the intersection of intervals for each trials with stimON
    stimTrials = [];
    for i = 1:numel(trials)
        t = InIntervals(behavior.stimON, trials{i}.timestamps);
        if sum(t) > 5 % I'm not sure why there are a few stim trials in the wrong direction
            stimTrials(i, 1) = 1;
            behavior.trialIDname{i, 2} = 'stimON';
            %trials{i}.manipulation = 'ON';
        else
            stimTrials(i, 1) = 0;
            behavior.trialIDname{i, 2} = 'stimOFF';
            %trials{i}.manipulation = 'OFF';
        end
        clear t;
    end
    for i = 1:numel(behavior.trialID)
        if behavior.trialID(i) == 1
            behavior.trialID(i, 2) = stimTrials(1);
        elseif behavior.trialID(i) == 2
            behavior.trialID(i, 2) = stimTrials(2);
        end
    end

end

%% Plots to check results
% needs improvement
if show_fig
    figure;
    plot(behavior.timestamps, behavior.position.linearized-min(behavior.position.linearized), 'k', 'LineWidth', 2);
    hold on;
    plot(behavior.timestamps, behavior.speed, 'r', 'LineWidth', 2);
    hold on;
    PlotIntervals(trials{1}.timestampsRun, 'color', 'b', 'alpha', .5);
    hold on;
    PlotIntervals(trials{2}.timestampsRun, 'color', 'g', 'alpha', .5);
    hold on;
    if ~isempty(manipulation) && exist([basepath, filesep, [basename, '.pulses.events.mat']], 'file')
        PlotIntervals(behavior.stimON, 'color', 'k', 'alpha', .5);
        hold on;
    end
    ylim([-1, max(behavior.position.linearized) - min(behavior.position.linearized)])
    saveas(gcf, [basepath, filesep, [basename, '.linearTrackBehavior.fig']]);
end

if remove_extra_fields
    behavior = rmfield(behavior, {'positionTrials', 'run', 'positionTrialsRun', 'speedTh'});
end

%% Generate output variables
if savemat
    if just_save_animal_behavior
        save([basepath, filesep, [basename, '.animal.behavior.mat']], 'behavior');
    else
        save([basepath, filesep, [basename, '.animal.behavior.mat']], 'behavior');
        save([basepath, filesep, [basename, '.trials.mat']], 'trials');
    end
end

end
