function speed_smooth = get_SmoothedSpeed(behavior,MergePoints,epochs,orderKalmanVel,do_smooth,varargin)

% USAGE
% speed_smooth = get_SmoothedSpeed(behavior,MergePoints,epoch,order)
% Calculates smoothed running speed using the Kalman Filter on
% recorded positions
%
%       calls KalmanVel.m
%
% INPUTS
%   behavior  - buzcode format behavior struct, to get positions
%   MergePoints - buzcode format MergePoints struct, to get the start and
%   end time of epochs
%   epochs - epoch numbers
%   orderKalmanVel - order of Kalman Velocity Filter (default 2)
%   do_smooth - gaussian smooth the positions? 1  = yes, 0 = no
%
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'	    - file path (default = pwd)
%     'smoothing_width'	- SD of gaussian used to smooth position (default = 5 for mice; recommend 10 for rats)
%     'saveMat'   		- logical (default: false) that updates file
%     'figopt'			- logical (default: false) that plots the raw and smoothed speed
%
% OUTPUT
%
%   speed_smooth - smoothed speed, aligned to behavior.timestamps
%
% Wenbo Tang, Jan 31, 2023
%% parse inputs
p=inputParser;
addParameter(p,'smoothing_width',5,@isnumeric);
addParameter(p,'figopt',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'basepath',pwd,@ischar);

parse(p,varargin{:});
smoothing_width = p.Results.smoothing_width;
figopt = p.Results.figopt;
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
%%
% Calculate speed based on 2D positions, instead of linearized, is
% recommended. As the animal may be actually moving, but not changing its
% linearized position. 

% get 2D position
positions = [behavior.timestamps',behavior.position.x',behavior.position.y'];
% linpos = [behavior.timestamps',behavior.position.linearized']; % you may use the linearized position, but not recommended
%%
% epoch loop
speed_smooth = [];% reset the matrix
for ep = epochs 
    % restrict to the current epoch
    epochtimes = MergePoints.timestamps(ep,:); % current epoch
   
    valid_id = find(positions(:,1) >= epochtimes(1) & positions(:,1) <= epochtimes(2)); % restrict to the current epoch
    positions_ep = positions(valid_id,:);
    
    %find valid position points
    validpos_id = find(~isnan(positions_ep(:,2)) & ~isnan(positions_ep(:,3)));
    positions_ep_val = positions_ep(validpos_id,:);

    %take care of the edges
    if positions_ep_val(end,1) < positions_ep(end,1) % the last pos value is missing
       positions_ep_val(end+1,:) = positions_ep(end,:);
    end

    if positions_ep_val(1,1) > positions_ep(1,1) % the first pos value is missing
       positions_ep_val = [positions_ep(1,:);positions_ep_val];
    end  

    % interpolation to fill the gaps
    posx_pad = interp1(positions_ep_val(:,1),positions_ep_val(:,2),positions_ep(:,1),'pchip');
    posy_pad = interp1(positions_ep_val(:,1),positions_ep_val(:,3),positions_ep(:,1),'pchip');
    
    if do_smooth
        % gaussian smooth the position signal
        posfilt = gaussian(smoothing_width, ceil(4*smoothing_width)); % gaussian smoothing for velocity filter
        posx_pad = filtfilt(posfilt,1,posx_pad);
        posy_pad = filtfilt(posfilt,1,posy_pad);
    end
    
    [~,~,~,vx,vy,~,~] = KalmanVel(posx_pad,posy_pad,positions_ep(:,1),orderKalmanVel);%calculate speed with Kalman filters
    speed_ep = sqrt(vx.^2+vy.^2);
    
    speed_smooth = [speed_smooth;speed_ep];
end
%%
if figopt % plot the result
    plot(behavior.timestamps, behavior.speed) % raw speed in file
    hold on
    plot(behavior.timestamps,speed_smooth','.')
end
