function good_idx = manual_trackerjumps(ts,x,y,StartofRec,EndofRec,basepath,session_paths, varargin)

% Manually cut out xy coordinates that are outside the bounds of your maze.
%
% These can be caused by unplugs or if the rat jumps out.
% If you do not remove these points, your ratemap will be messed up.
%
% Input:
%       ts
%       x
%       y
%       StartofRec: ts indicating the start points of your event
%       EndofRec: ts indicating the end points of your event
%       basepath: path to session folder
%       session_paths: array of full paths for each subsession (for video overlay)
%
% Optional Parameters:
%       'video_overlay': logical, overlay tracking on video frame (default: false)
%       'darkmode': logical, use dark theme for plots (default: true)
%       'save_boundary_file': logical, save boundary selection (default: false)
%
% Copyright (C) 2018 Ryan E Harvey

% 1 for elminating the points outside drawn shape & 0 for inside
p = inputParser;
addParameter(p,'darkmode',true,@islogical);
addParameter(p,'restic_dir',1,@isnumeric);
addParameter(p,'axis_equal',true,@islogical);
addParameter(p,'alpha',.2,@isnumeric);
addParameter(p,'add_scatter',true,@islogical);
addParameter(p,'save_boundary_file',false,@islogical);
addParameter(p,'video_overlay',false,@islogical);

parse(p,varargin{:});
darkmode = p.Results.darkmode;
restic_dir = p.Results.restic_dir;
axis_equal = p.Results.axis_equal;
alpha = p.Results.alpha;
add_scatter = p.Results.add_scatter;
save_boundary_file = p.Results.save_boundary_file;
video_overlay = p.Results.video_overlay;
video_frame_time = [];

savets=[];
for i=1:length(StartofRec)
    % index out each event
    xtemp=x(ts>=StartofRec(i) & ts<=EndofRec(i));
    ytemp=y(ts>=StartofRec(i) & ts<=EndofRec(i));
    tstemp=ts(ts>=StartofRec(i) & ts<=EndofRec(i));
    
    % Check if we have valid data
    if isempty(xtemp) || isempty(ytemp)
        warning('Event %d has no data points', i);
        continue;
    end
    
    % Remove NaN values for visualization
    valid_idx = ~isnan(xtemp) & ~isnan(ytemp);
    if sum(valid_idx) == 0
        warning('Event %d has only NaN values', i);
        continue;
    end
    
    fprintf('Event %d: %d total points, %d valid points\n', i, length(xtemp), sum(valid_idx));

    % Auto-detect video path if video_overlay is enabled but path is empty
    current_video_path = '';
    basename = basenameFromBasepath(session_paths{i});
    if p.Results.video_overlay && isempty(current_video_path)
        fprintf('Video overlay enabled, searching for video file...\n');
        
        % Look for video files in MergePoints folders first
        if exist(fullfile(session_paths{i}, [basename, '.MergePoints.events.mat']), 'file')
            load(fullfile(session_paths{i}, [basename, '.MergePoints.events.mat']), 'MergePoints');
            
            % Check subsession folders for MP4/AVI files
            for k = 1:length(MergePoints.foldernames)
                mp4_files = dir(fullfile(session_paths{i}, MergePoints.foldernames{k}, '*.mp4'));
                if ~isempty(mp4_files)
                    current_video_path = fullfile(session_paths{i}, MergePoints.foldernames{k}, mp4_files(1).name);
                    fprintf('Found video: %s\n', current_video_path);
                    break;
                end
                
                avi_files = dir(fullfile(session_paths{i}, MergePoints.foldernames{k}, '*.avi'));
                if ~isempty(avi_files)
                    current_video_path = fullfile(session_paths{i}, MergePoints.foldernames{k}, avi_files(1).name);
                    fprintf('Found video: %s\n', current_video_path);
                    break;
                end
            end
        end
        
        % If still not found, check main session folder
        if isempty(current_video_path)
            mp4_files = dir(fullfile(session_paths{i}, '*.mp4'));
            if ~isempty(mp4_files)
                current_video_path = fullfile(session_paths{i}, mp4_files(1).name);
                fprintf('Found video: %s\n', current_video_path);
            else
                avi_files = dir(fullfile(session_paths{i}, '*.avi'));
                if ~isempty(avi_files)
                    current_video_path = fullfile(session_paths{i}, avi_files(1).name);
                    fprintf('Found video: %s\n', current_video_path);
                end
            end
        end
        
        if isempty(current_video_path)
            fprintf('Warning: No video file found for overlay. Proceeding without video.\n');
        end
    end
    
    % use the gui to cut out points
    [~,~,in]=restrictMovement(xtemp,ytemp,restic_dir,darkmode,axis_equal,...
        alpha,add_scatter,video_overlay,current_video_path,video_frame_time,tstemp);
    % save the ts where the tracker error exists
    savets=[savets,tstemp(in)];
end

% locate the index for each tracker error
good_idx=ismember(ts,savets);

% save that index to your session folder so you won't have to do this again
% each time you run your data
if save_boundary_file
    basename = basenameFromBasepath(basepath);
    save(fullfile(basepath,[basename,'.restrictxy.mat']),'good_idx')
end
end

function [x,y,in]=restrictMovement(x,y,direction,darkmode,axis_equal,alpha,add_scatter,video_overlay,video_path,video_frame_time,timestamps)
% restrictMovement allows you to draw a line around xy coordinates in order
% to eliminate certain points you don't want...ie when the rat jumps out of
% maze or tracker errors.
%
% Also, I have included a direction input argument so you have either
% restrict outside or inside points. This is valuable if you have a maze
% like a circular track where the rat could jump out away or towards the
% center of the maze.
%
% Input         x,y: coordinates
%         direction: 1 (default) to remove outside points; 0 to remove inside points
%         video_overlay: logical to overlay video frame as background
%         video_path: path to video file
%         video_frame_time: specific time to extract frame, or empty for middle frame
%         timestamps: timestamps for current event
%
%
% Output        x,y: retained coordinates
%                in: logical of which coordinates were retained (so you can index ts etc.)
%
%
% Ryan Harvey

% check inputs
if nargin<3
    direction=1;
end
if nargin<8
    video_overlay=false;
end
if nargin<9
    video_path='';
end
if nargin<10
    video_frame_time=[];
end
if nargin<11
    timestamps=[];
end

% Validate input data
if isempty(x) || isempty(y)
    error('Empty x or y coordinates');
end

% Remove NaN values for visualization
valid_idx = ~isnan(x) & ~isnan(y);
if sum(valid_idx) == 0
    error('All coordinates are NaN');
end

x_plot = x(valid_idx);
y_plot = y(valid_idx);

fprintf('Plotting %d valid points (out of %d total)\n', sum(valid_idx), length(x));
fprintf('X range: [%.2f, %.2f], Y range: [%.2f, %.2f]\n', ...
    min(x_plot), max(x_plot), min(y_plot), max(y_plot));

% Load video frame if requested
video_frame = [];
if video_overlay && ~isempty(video_path) && exist(video_path, 'file')
    try
        fprintf('Loading video frame from: %s\n', video_path);
        
        % Determine which frame to extract
        if isempty(video_frame_time) && ~isempty(timestamps)
            % Use middle timestamp of current event
            video_frame_time = median(timestamps) - timestamps(1);
            fprintf('Using median timestamp: %.3f s\n', video_frame_time);
        elseif isempty(video_frame_time)
            % Use middle of video
            video_frame_time = [];
            fprintf('Using middle frame of video\n');
        else
            fprintf('Using specified timestamp: %.3f s\n', video_frame_time);
        end
        
        % Try to read video frame
        try
            v = VideoReader(video_path);
            
            % Check if video_frame_time exceeds video duration
            if ~isempty(video_frame_time) && video_frame_time > v.Duration
                fprintf('Warning: Requested time (%.3f s) exceeds video duration (%.3f s). Skipping video overlay.\n', ...
                    video_frame_time, v.Duration);
                % use previous frame if available
                if exist('video_frame', 'var') && ~isempty(video_frame)
                    fprintf('Using previously loaded frame\n');
                else
                    video_frame = [];
                end
            elseif ~isempty(video_frame_time)
                v.CurrentTime = video_frame_time;
                if hasFrame(v)
                    video_frame = readFrame(v);
                    fprintf('Successfully loaded video frame\n');
                else
                    fprintf('Warning: No frame available at specified time\n');
                    % use previous frame if available
                    if exist('video_frame', 'var') && ~isempty(video_frame)
                        fprintf('Using previously loaded frame\n');
                    else
                        video_frame = [];
                    end
                end
            else
                v.CurrentTime = v.Duration / 2; % Middle of video
                if hasFrame(v)
                    video_frame = readFrame(v);
                    fprintf('Successfully loaded video frame\n');
                else
                    fprintf('Warning: No frame available at specified time\n');
                    % use previous frame if available
                    if exist('video_frame', 'var') && ~isempty(video_frame)
                        fprintf('Using previously loaded frame\n');
                    else
                        video_frame = [];
                    end
                end
            end
        catch video_err
            fprintf('Warning: VideoReader failed: %s\n', video_err.message);
            fprintf('Attempting ffmpeg method...\n');
            
            % Get video duration using ffprobe for validation
            if ~isempty(video_frame_time)
                duration_cmd = sprintf('ffprobe -v error -show_entries format=duration -of default=noprint_wrappers=1:nokey=1 "%s" 2>/dev/null', video_path);
                [status_dur, duration_str] = system(duration_cmd);
                if status_dur == 0
                    video_duration = str2double(strtrim(duration_str));
                    if video_frame_time > video_duration
                        fprintf('Warning: Requested time (%.3f s) exceeds video duration (%.3f s). Skipping video overlay.\n', ...
                            video_frame_time, video_duration);
                        video_frame = [];
                        % Skip ffmpeg extraction
                        continue_to_ffmpeg = false;
                    else
                        continue_to_ffmpeg = true;
                    end
                else
                    % Could not get duration, try ffmpeg anyway
                    continue_to_ffmpeg = true;
                end
            else
                continue_to_ffmpeg = true;
            end
            
            % Fallback: try to extract frame using system ffmpeg
            if continue_to_ffmpeg
                temp_frame_path = tempname;
                temp_frame_path = [temp_frame_path, '.jpg'];
                
                if ~isempty(video_frame_time)
                    cmd = sprintf('ffmpeg -i "%s" -ss %.3f -vframes 1 -y "%s" 2>/dev/null', ...
                        video_path, video_frame_time, temp_frame_path);
                else
                    cmd = sprintf('ffmpeg -i "%s" -vframes 1 -y "%s" 2>/dev/null', ...
                        video_path, temp_frame_path);
                end
                
                [status, ~] = system(cmd);
                if status == 0 && exist(temp_frame_path, 'file')
                    video_frame = imread(temp_frame_path);
                    delete(temp_frame_path);
                    fprintf('Successfully extracted frame using ffmpeg\n');
                else
                    fprintf('Warning: ffmpeg extraction also failed\n');
                    if exist(temp_frame_path, 'file')
                        delete(temp_frame_path);
                    end
                end
            end
        end
        
    catch err
        fprintf('Warning: Could not load video frame: %s\n', err.message);
        video_frame = [];
    end
elseif video_overlay && isempty(video_path)
    fprintf('Warning: video_overlay requested but no video_path provided\n');
elseif video_overlay && ~exist(video_path, 'file')
    fprintf('Warning: video file not found: %s\n', video_path);
end

% set up figure
figure;

% Display video frame as background if available
if ~isempty(video_frame)
    % Display video frame
    imshow(video_frame);
    hold on;
    
    % Get image dimensions for coordinate scaling
    [img_height, img_width, ~] = size(video_frame);
    
    % Scale tracking coordinates to match image dimensions if needed
    % Assume tracking coordinates are already in pixel coordinates
    x_scaled = x_plot;
    y_scaled = y_plot;
    
    % Plot tracking data over video
    if darkmode
        plot(x_scaled, y_scaled, 'Color', [1,1,1,alpha], 'LineWidth', 1.5);
        if add_scatter
            scatter(x_scaled, y_scaled, 8, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    else
        plot(x_scaled, y_scaled, 'Color', [1,0,0,alpha], 'LineWidth', 1.5);
        if add_scatter
            scatter(x_scaled, y_scaled, 8, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    end
    
    title('Click around the points you want to keep (Video overlay)')
    
    % Set axis limits to match image
    xlim([0, img_width]);
    ylim([0, img_height]);
    axis equal;
    
else
    % Original plotting without video overlay
    if darkmode
        plot(x_plot,y_plot,'Color',[1,1,1,alpha]);hold on
        if add_scatter
            scatter(x_plot,y_plot,3,'w','filled');hold on
        end
        title('Click around the points you want to keep')
        xlabel('X')
        ylabel('Y')
        axis tight
        if axis_equal
            axis equal
        end
        darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
    else
        plot(x_plot,y_plot,'Color',[0,0,0,alpha]);hold on
        if add_scatter
            scatter(x_plot,y_plot,3,'k','filled');hold on
        end
        title('Click around the points you want to keep')
        xlabel('X')
        ylabel('Y')
        axis tight
        if axis_equal
            axis equal
        end
    end
end
disp('PRESS "ENTER" TO EXIT')
i=1;

% let the user click around the coordinates
while true
    [X,Y]=ginput(1);
    if isempty(X)
        break
    end
    corners(i,:)=[X,Y];
    plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
    i=i+1;
end

% remove points outside or inside the shape
in=inpolygon(x,y,corners(:,1),corners(:,2));
if direction==1
    x=x(in);
    y=y(in);
else
    x=x(~in);
    y=y(~in);
    in=~in;
end
close
end