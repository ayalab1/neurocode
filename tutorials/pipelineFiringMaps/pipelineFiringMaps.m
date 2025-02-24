%% Basic pipeline for calculating rate maps
%---------------------------------------------------------------%
%  This is the main script for a step-by-step tutorial for      %
%  calculating rate maps for each cells                         % 
%  -- Wenbo Tang (Jan 31, 2023)                                 %
%---------------------------------------------------------------%
clc
clear all
close all
%%
% information of your data to be analyzed
filepath = '/Volumes/AYALab_Data/Sample_data/';
animalprefix = 'HYC3';
day = 26;
daystring = num2str(day);
basepath = [filepath,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
epochs = [1,3]; % behavioral epochs, not sleep epochs
conditions = 2; % number of trajectory types/ conditions (e.g., 2 for a typical linear track, inbound and outbound)
%%
%----------set parameters-----------%
CellType = "Pyramidal Cell"; % cell types to be analyzed; set to empty to load all
brainRegion = "CA1"; % brain regions to be analyzed; set to empty to load all
speedthresh = 2; % apply speed threshold (in cm/s); default threshold = 2 for mice, 5 for rats; set to empty if no threshold used
savedata = true; % save results?

% select the analysis to apply
recalculate_speed = true; % recalculate speed using Kalman filters
use_smoothed_speed = true; % use the smooth_speed, instead of raw speed
% **** Note: place-field statistics (calculate_ratemap_stats) will only perform when calculate_ratemap is true
calculate_2D_ratemap = true; % calculate 2D rate maps
calculate_2D_ratemap_stats = true; % calculate statistic of 2D place fields
calculate_linear_ratemap = false; % calculate linearized rate maps
calculate_linear_ratemap_stats = false; % calculate statistic of linearized place fields
%%
%----------set parameters for analysis-----------%
% only if recalculate_speed = 1, define the parameters below; otherwise, skip
if recalculate_speed
    do_smooth = 1; % gaussian smooth the position
    smoothing_width = 5; % SD of gaussian used to smooth position, recommend 5 for mice, 10 for rats
    orderKalmanVel = 2; % order of the Kalman filters
    doPlot = true; % plot the smoothed speed
end

if calculate_2D_ratemap || calculate_linear_ratemap
    % parameters for rate map calculation
    smooth = 2; % smoothing size in bins (0 = no smoothing, default = 2)
    nBins = 100; % number of bins for linearized rate mpas (default = [50 50])
    nBins_2D = 150; % number of horizontal and vertical bins for 2D rate maps (default = [50 50])
    maxGap = 0.1; % z values recorded during time gaps between successive (x,y)
    minTime = 0.02; % minimum time spent in each bin (in s, default = 0)
    mode = 'discard'; % 'interpolate' to interpolate missing points (< minTime)
    maxDistance = 5; % maximal distance for interpolation (default = 5)
    % Optional:
    % 'type'            three letters (one for X, one for Y and one for Z) indi-
    %                   cating which coordinates are linear ('l') and which are
    %                   circular ('c') - for 1D data, only two letters are used
    %                   (default 'lll')
end

if calculate_2D_ratemap_stats || calculate_linear_ratemap_stats
    smoothing_width = 5; % SD of gaussian used to smooth position, recommend 5 for mice, 10 for rats
    nShuffles = 100; % number of shuffles to calculate the significance of spatial information
    threshold = 0.15; % values above threshold*peak belong to the field (default = 0.15)
    minSize = 0.05; % fields smaller than this percentage of the maze size are considered spurious and ignored
    maxSize = 0.50; % fields larger than this percentage of the maze size are considered noise and ignored
    sepEdge = 0; % fields with maximum Firing Rate closer to the edges less than this percentage of the maze size are ignored
    minPeak = 1; % peaks smaller than this size are considered spurious and ignored (default = 1 Hz)
    minPeak2nd = 0.60;%  for secondary place fields, peaks smaller than this  percentage of maximum Firing Rate along the maze are
% considered spurious and ignored (default 0.60)
    doPlot = false; % Make plots (default: true) 
    posTrials_smooth = true; % smooth out the jumping points of linearized positions from tracking errors
    if posTrials_smooth
        gapdis  = 10;% the threshold for distance between tracking errors and actual position, in cm 
    end
end

%% load files for preprocessing
load(fullfile(basepath,[prefix,'.animal.behavior.mat'])) % behavior info to get position, linearized position, speed and etc
load(fullfile(basepath,[prefix,'.MergePoints.events.mat'])) % Session info to get start and end time of each epoch
spikes = importSpikes('basepath',basepath,'CellType',CellType,'brainRegion',brainRegion); % load spikes
load(fullfile(basepath,[prefix, '.cell_metrics.cellinfo.mat'])); % cell information
%% recalculate speed?
if recalculate_speed
    speed_smooth = get_SmoothedSpeed(behavior,MergePoints,epochs,orderKalmanVel,do_smooth,...
        'smoothing_width',smoothing_width,'figopt', doPlot);
    behavior.speed_smooth = speed_smooth';
    if savedata
       save(fullfile(basepath,[prefix,'.animal.behavior.mat']),'behavior') % update behavior info
    end
end
% use the smoothed speed?
if use_smoothed_speed
    if ~isfield(behavior,'speed_smooth')
        warning('No "speed_smooth" calculated in the behavior file, use raw speed instead')
        vel  = behavior.speed'; % raw running speed
    elseif isempty(behavior.trials)
        vel  = behavior.speed_smooth'; % smoothed running speed (recommended)
    end
else
    vel  = behavior.speed'; % raw running speed
end
%% calculate 2D rate maps
if calculate_2D_ratemap
    % gather speed-thresholded positions
    positions = [behavior.timestamps',behavior.position.x',behavior.position.y',vel]; %[t,x,y,vel]
    
    % calculate 2D rate maps for each epoch; combining epoches is not an
    % option for now
    cellsample_count = 0; % count the number of cell samples (unit in each condition)
    for epoch = epochs
        % restrict to the current epoch
        epochtimes = MergePoints.timestamps(epoch,:); % current epoch

        valid_id = find(positions(:,1) >= epochtimes(1) & positions(:,1) <= epochtimes(2)); % restrict to the current epoch
        positions_ep = positions(valid_id,:);

        % apply speed threshold
        if ~isempty(speedthresh)
            validtime_id = find(positions_ep(:,4) >= speedthresh);
            temp = positions_ep(validtime_id,:);
            positions_valid = temp;
        else
            positions_valid = positions_ep;
        end
        %%
        %-------calculate 2D rate maps for each cells--------%
        cellnum = length(spikes.times); % number of neurons
        for unit = 1:length(spikes.times)
            
           
            temp = Map(positions_valid(:,1:3),spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                'nBins',nBins_2D,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
            
            if ~isempty(temp.z)
                cellsample_count = cellsample_count +1;


                % set low occ bins to -1
                invalidid = find(temp.time < minTime);
                temp.z(invalidid) = -1;


                %--------do place field detection and statistics?-------%
                if calculate_2D_ratemap_stats
                    % for 2D rate maps, field detetction is not implemented
                        validid = find(temp.time >= minTime);
                        temp_valid.z = temp.z(validid);
                        temp_valid.count = temp.count(validid);
                        temp_valid.time = temp.time(validid);

                        % set low occ bins to -1
                        stats = MapStats2D_basic(temp_valid);
                        % shuffles to get significant for specificity (or
                        % spatial information)
                        specificity_shuf = nan(nShuffles,1);
                        maxTime = max(positions_valid(:,1));
                        for i = 1:nShuffles
                            shifted = sortrows([rem(positions_valid(:,1)+rand(1)*maxTime,maxTime) positions_valid(:,2:3)]);
                            temp_shuffled = Map(shifted,spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                    'nBins',nBins_2D,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);

                            temp_valid_shuf.z = temp_shuffled.z(validid);
                            temp_valid_shuf.count = temp_shuffled.count(validid);
                            temp_valid_shuf.time = temp_shuffled.time(validid);
                            stats_shuf = MapStats2D_basic(temp_valid_shuf);
                            specificity_shuf(i) = stats_shuf.specificity;
                        end
                        specificity_pval = 1-sum(specificity_shuf < stats.specificity)./sum(~isnan(specificity_shuf));
                end


                %%
                % organizing into struct
                % rate maps
                firingMaps_2D.ratemaps{cellsample_count} = temp.z;
                firingMaps_2D.countMaps{cellsample_count} = temp.count;
                firingMaps_2D.occupancy{cellsample_count} = temp.time;

                % place field statistics
                if calculate_2D_ratemap_stats
                    firingMaps_2D.info.specificity(cellsample_count) = stats.specificity;
                    firingMaps_2D.info.specificity_pval(cellsample_count) = specificity_pval;
                    firingMaps_2D.info.informationPerSpike(cellsample_count) = stats.informationPerSpike;
                    firingMaps_2D.info.informationPerSec(cellsample_count) = stats.informationPerSec;
                    firingMaps_2D.info.sparsity(cellsample_count) = stats.sparsity;
                    firingMaps_2D.info.selectivity(cellsample_count) = stats.selectivity;
                end

                % analysis parameters
                firingMaps_2D.info.params_smooth(cellsample_count) = smooth;
                firingMaps_2D.info.params_minTime(cellsample_count)  = minTime;
                firingMaps_2D.info.params_nBins(cellsample_count)  = nBins_2D;
                firingMaps_2D.info.params_x{cellsample_count}  = temp.x;
                firingMaps_2D.info.params_y{cellsample_count}  = temp.y;
                firingMaps_2D.info.params_maxGap(cellsample_count)  = maxGap;
                firingMaps_2D.info.params_mode{cellsample_count}  = mode;
                firingMaps_2D.info.params_maxDistance(cellsample_count)  = maxDistance;

                % data info
                firingMaps_2D.info.UID(cellsample_count) = spikes.UID(unit);
                cellid = find(cell_metrics.UID == spikes.UID(unit));
                firingMaps_2D.info.brainRegion{cellsample_count} = cell_metrics.brainRegion{cellid};
                firingMaps_2D.info.putativeCellType{cellsample_count} = cell_metrics.putativeCellType{cellid};
                firingMaps_2D.info.epochname{cellsample_count} = behavior.epochs{epoch}.name;
                firingMaps_2D.info.epoch_num(cellsample_count) = epoch -1; %start with 0
                firingMaps_2D.info.epoch_duration(cellsample_count) = behavior.epochs{epoch}.stopTime - behavior.epochs{epoch}.startTime;
                firingMaps_2D.info.epoch_start(cellsample_count) = behavior.epochs{epoch}.startTime;
                firingMaps_2D.info.epoch_stop(cellsample_count) = behavior.epochs{epoch}.stopTime;
            end
        end
    end
    %--------save the 2D place fields-------%
    if savedata && cellsample_count > 0
       save(fullfile(basepath, [basenameFromBasepath(basepath),'.firingMapsAvg_2D.cellinfo.mat']),'firingMaps_2D');
    end
end

%% calculate 1D linearized rate maps
% **** Note: you need to log the state transition of each trajectory type based
% on the layout of your linearizaion.
if calculate_linear_ratemap
    %------get to know your linearization config------%
    % **** Note: this section needs to be tailored for each behavior task. 
    % info of linearization segments
    linpos_states = unique(behavior.states(~isnan(behavior.states)));% number of segments/ states; state starts from 0
    state_distance = [];% length of each segment
    state_posxy = [];% [x,y] coordinates of the start and end positions of each segment
    figure,% plot a linearized skeleton to get an idea of the general layout 
    for s = linpos_states
        stateid = find(behavior.states == s);
        [startdis,startid] = min(behavior.position.linearized(stateid));
        startid = stateid(startid);
        [enddis,endid] = max(behavior.position.linearized(stateid));
        endid = stateid(endid);
        start_posxy = [behavior.position.x(startid),behavior.position.y(startid)];
        end_posxy = [behavior.position.x(endid),behavior.position.y(endid)];

        state_distance = [state_distance;startdis,enddis];
        state_posxy = [state_posxy;start_posxy,end_posxy];  

        plot([start_posxy(1),end_posxy(1)],[start_posxy(2),end_posxy(2)])
        hold on
        text(start_posxy(1) +(end_posxy(1)-start_posxy(1))/2,...
            start_posxy(2) +(end_posxy(2)-start_posxy(2))/2,['\leftarrow ',num2str(s)]);
    end
    hold off
    
    %% please manually log the state transition of each trajectory type (or conditions) here
    % demo for a standard T-maze (2 conditions, left-side vs right-side trajectory)
    keyboard; % a reminder for checking the state transition matrix
    % state transition for the 2 trajectory types
    statematrix{1} = [0,4,7,8,2];%left-side
    statematrix{2} = [0,3,5,6,1];%right-side

    centerseg_incm = state_distance(1,2);%choice point
    rightend_incm = state_distance(2,2);%end point of right trajectories
    leftend_incm = state_distance(3,2);%end point of left trajectories
    leftend_incm = leftend_incm-rightend_incm + centerseg_incm;
    %%
    %-----get linearized position for each trial/ trajectory----%
    % **** Note: this section needs to be tailored for each behavior task. 
    % we need the start and end time for each trial. Use the position
    % pipeline the define trials first.
    if ~isfield(behavior,'trials')
        error('No "trials" defined in the behavior file.')
    elseif isempty(behavior.trials)
        error('No "trials" defined in the behavior file.')
    end
    % we need the condition/ trajectory type for each trial. Use the position
    % pipeline the define trial conditions first.
    if ~isfield(behavior,'trialConds')
        error('No "trial conditions" defined in the behavior file.')
    elseif isempty(behavior.trialConds)
        error('No "trial conditions" defined in the behavior file.')
    end
    
    % get linearized position for each trials
    for iCond = 1:conditions
        posTrials{iCond} = [];
    end
    
    trialnum = length(behavior.trials(:,1));
    for tr = 1:trialnum
        currenttrial_time = behavior.trials(tr,:);
        trialid = find(behavior.timestamps >= currenttrial_time(1) & behavior.timestamps <= currenttrial_time(2));
        current_posTrial = [behavior.timestamps(trialid)',behavior.position.linearized(trialid)']; % raw posTrials
        if behavior.trialConds(tr) == 1 % left-side trial
           % correct projection error around choice point
           errorid = find(current_posTrial(:,2) > centerseg_incm & current_posTrial(:,2) <= rightend_incm);
           current_posTrial(errorid,2)  =  current_posTrial(errorid,2) + rightend_incm;

           % convert to the distance to start
           sideseg_id = find(current_posTrial(:,2) > rightend_incm);
           current_posTrial(sideseg_id,2) = current_posTrial(sideseg_id,2) - rightend_incm + centerseg_incm; % corrected posTrials

           if posTrials_smooth %smooth out the jumping points from tracking errors
               validid = find(~isnan(current_posTrial(:,2)));
               posfilt = gaussian(smoothing_width, ceil(4*smoothing_width)); % gaussian smoothing for velocity filter
               linpos_smooth = filtfilt(posfilt,1,current_posTrial(validid,2));
               linpos_diff = abs(linpos_smooth - current_posTrial(validid,2));
               invalid_id  = find(linpos_diff > gapdis);
               current_posTrial(validid(invalid_id),2) = linpos_smooth(invalid_id);
           end

           % save the result
           posTrials{1} = [posTrials{1};current_posTrial];
        elseif behavior.trialConds(tr) == 2
           % correct projection error around choice point
           errorid = find(current_posTrial(:,2) > rightend_incm);
           current_posTrial(errorid,2)  =  current_posTrial(errorid,2) - rightend_incm;
           
           if posTrials_smooth %smooth out the jumping points from tracking errors
               validid = find(~isnan(current_posTrial(:,2)));
               posfilt = gaussian(smoothing_width, ceil(4*smoothing_width)); % gaussian smoothing for velocity filter
               linpos_smooth = filtfilt(posfilt,1,current_posTrial(validid,2));
               linpos_diff = abs(linpos_smooth - current_posTrial(validid,2));
               invalid_id  = find(linpos_diff > gapdis);
               current_posTrial(validid(invalid_id),2) = linpos_smooth(invalid_id);
           end
           % save the result
           posTrials{2} = [posTrials{2};current_posTrial];
        end
    end

    % apply speed threshold
    if ~isempty(speedthresh)
        validtime = behavior.timestamps(find(vel >= speedthresh)); % find running time periods
        for iCond = 1:conditions
            [~,~,linpos_id] = intersect(validtime,posTrials{iCond}(:,1));
            posTrials{iCond} = posTrials{iCond}(linpos_id,:);
        end
    end
    %%
    %-------calculate linearized rate maps and their statistics for each cells--------%
    cellsample_count = 0; % count the number of cell samples (unit in each condition)
    for epoch = epochs
        % restrict to the current epoch
        epochtimes = MergePoints.timestamps(epoch,:); % current epoch
        for iCond = 1:conditions
            valid_id = find(posTrials{iCond}(:,1) >= epochtimes(1) & posTrials{iCond}(:,1) <= epochtimes(2)); % restrict to the current epoch
            posTrials_ep{iCond} = posTrials{iCond}(valid_id,:);
        end
       
        %%
        %-------calculate linearized rate maps for each cells--------%
        cellnum = length(spikes.times); % number of neurons
        for unit = 1:length(spikes.times)
            % different conditions
            for c = 1:conditions

                temp = Map(posTrials_ep{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                    'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                if ~isempty(temp.z)
                    cellsample_count = cellsample_count +1;
                    %--------do place field detection and statistics?-------%
                    if calculate_linear_ratemap_stats


                        % Default values
                        stats.x = NaN;
                        stats.field = [];
                        stats.size = 0;
                        stats.peak = 0;
                        stats.mean = 0;
                        stats.fieldX = [NaN NaN];

                        stats.specificity = nan;
                        stats.informationPerSpike = nan;
                        stats.informationPerSec = nan;
                        stats.sparsity = nan;
                        stats.selectivity = nan;
                        specificity_pval = nan;

                    
                        stats = MapStats1D(temp,'threshold',threshold,'minSize',minSize,...
                    'maxSize',maxSize,'minPeak',minPeak,'minPeak2nd',minPeak2nd,'sepEdge',sepEdge,...
                    'doPlot',doPlot);
                
                        % shuffles to get significant for specificity (or
                        % spatial information)
                        specificity_shuf = nan(nShuffles,1);
                        maxTime = max(posTrials_ep{c}(:,1));
                        for i = 1:nShuffles
                            shifted = sortrows([rem(posTrials_ep{c}(:,1)+rand(1)*maxTime,maxTime) posTrials_ep{c}(:,2)]);
                            temp_shuffled = Map(shifted,spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                    'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                
                           
                        stats_shuf = MapStats1D(temp_shuffled,'threshold',threshold,'minSize',minSize,...
                                            'maxSize',maxSize,'minPeak',minPeak,'minPeak2nd',minPeak2nd,'sepEdge',sepEdge,...
                                            'doPlot',false);
                            specificity_shuf(i) = stats_shuf.specificity;
                        end
                        specificity_pval = 1-sum(specificity_shuf < stats.specificity)./sum(~isnan(specificity_shuf));
                    end

                    % set low occ bins to -1
                    invalidid = find(temp.time < minTime);
                    temp.z(invalidid) = -1;
                    %%
                    % organizing into struct
                    % rate maps
                    firingMaps_linear.ratemaps{cellsample_count} = temp.z;
                    firingMaps_linear.countMaps{cellsample_count} = temp.count;
                    firingMaps_linear.occupancy{cellsample_count} = temp.time;

                    % place field statistics
                    if calculate_linear_ratemap_stats
                        firingMaps_linear.info.specificity(cellsample_count) = stats.specificity;
                        firingMaps_linear.info.specificity_pval(cellsample_count)= specificity_pval;
                        firingMaps_linear.info.informationPerSpike(cellsample_count) = stats.informationPerSpike;
                        firingMaps_linear.info.informationPerSec(cellsample_count) = stats.informationPerSec;
                        firingMaps_linear.info.sparsity(cellsample_count) = stats.sparsity;
                        firingMaps_linear.info.selectivity(cellsample_count) = stats.selectivity;

                        firingMaps_linear.info.field_x{cellsample_count} = stats.x;
                        firingMaps_linear.info.field{cellsample_count} = stats.field;
                        firingMaps_linear.info.field_size{cellsample_count} = stats.size;
                        firingMaps_linear.info.field_peakrate{cellsample_count} = stats.peak;
                        firingMaps_linear.info.field_fieldX{cellsample_count} = stats.fieldX;
                    end

                    % analysis parameters
                    firingMaps_linear.info.params_smooth(cellsample_count) = smooth;
                    firingMaps_linear.info.params_minTime(cellsample_count)  = minTime;
                    firingMaps_linear.info.params_nBins(cellsample_count)  = nBins;
                    firingMaps_linear.info.params_x{cellsample_count}  = temp.x;
                    firingMaps_linear.info.params_y{cellsample_count}  = temp.y;
                    firingMaps_linear.info.params_maxGap(cellsample_count)  = maxGap;
                    firingMaps_linear.info.params_mode{cellsample_count}  = mode;
                    firingMaps_linear.info.params_maxDistance(cellsample_count)  = maxDistance;

                    firingMaps_linear.info.params_threshold(cellsample_count) = threshold;
                    firingMaps_linear.info.params_minSize(cellsample_count)  = minSize;
                    firingMaps_linear.info.params_maxSize(cellsample_count)  = maxSize;
                    firingMaps_linear.info.params_sepEdge(cellsample_count)  = sepEdge;
                    firingMaps_linear.info.params_minPeak(cellsample_count)  = minPeak;
                    firingMaps_linear.info.params_minPeak2nd(cellsample_count)  = minPeak2nd;

                    % data info
                    firingMaps_linear.info.UID(cellsample_count) = spikes.UID(unit);
                    cellid = find(cell_metrics.UID == spikes.UID(unit));
                    firingMaps_linear.info.brainRegion{cellsample_count} = cell_metrics.brainRegion{cellid};
                    firingMaps_linear.info.putativeCellType{cellsample_count} = cell_metrics.putativeCellType{cellid};
                    firingMaps_linear.info.epochname{cellsample_count} = behavior.epochs{epoch}.name;
                    firingMaps_linear.info.epoch_num(cellsample_count) = epoch -1; %start with 0
                    firingMaps_linear.info.epoch_duration(cellsample_count) = behavior.epochs{epoch}.stopTime - behavior.epochs{epoch}.startTime;
                    firingMaps_linear.info.epoch_start(cellsample_count) = behavior.epochs{epoch}.startTime;
                    firingMaps_linear.info.epoch_stop(cellsample_count) = behavior.epochs{epoch}.stopTime;
                end
            end
        end
    end
    %--------save the linearized place fields-------%
    if savedata && cellsample_count > 0
       save(fullfile(basepath, [basenameFromBasepath(basepath),'.firingMapsAvg_linear.cellinfo.mat']),'firingMaps_linear');
    end
end

    
    
   
   
    



