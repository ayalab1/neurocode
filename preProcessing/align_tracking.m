function delay = align_tracking(timestamps, positions, events, varargin)
% Heuristic to align behavior tracking trajectory with detected pulses
% (e.g. from IR barriers, etc.), assuming that events always happen in a
% fixed location. The function returns the delay of the trajectory relative
% to the events that minimizes the jitter of event locations. (This
% implicitly assumes that frames are dropped in the beginning or end of the
% recording and that timestamps between frames are linear in between).
% 
% OUTPUTS:
% delay     double  Delay (in seconds) of the trajectory timestamps
%                   relative to ephys evt timestamps so that
%                   behavior.timestamsps + delay = correct timestamps.

%% Input handling
p = inputParser;

addParameter(p,'fs',30,@isnumeric); % Behavior acquisition rate
addParameter(p,'showresults',true,@islogical);
addParameter(p,'maxrange',120,@islogical); % Max range of shift in seconds

parse(p,varargin{:})

fs = p.Results.fs;
showresults = p.Results.showresults;
maxrange = p.Results.maxrange;

%% Setup
shifts = round(-maxrange*fs):1:round(maxrange*fs);

%% Find behavior timestamps corresponding to events
for e = 1:length(events)
    % Reject evts that are outside of tracking timestamps.
    if min(abs(timestamps-events(e))) <= 1/fs
    	[~,evidcs(e)] = min(abs(timestamps-events(e)));
    else
        evidcs(e) = NaN;
    end
end
evidcs = evidcs(~isnan(evidcs));

for s = 1:length(shifts)
    idcs = evidcs - shifts(s); % Note: here evts are shifted rel to ts, unlike in the later application for correction.
    idcs = idcs(idcs>0 & idcs<length(timestamps));
    var(s) = std(positions(idcs));
end

% Get optimum delay
[~,delay] = min(var);
delay = shifts(delay)/fs;

%% Plotting
if showresults
    figure;
    
    subplot(1,3,1);
    hold on
    title('Original recording');
    plot(timestamps,positions);
    scatter(timestamps(evidcs),positions(evidcs));
    
    subplot(1,3,2);
    hold on
    title(sprintf('Optimal delay = %.2f s',delay));
    plot(shifts/fs,var);
    xlabel('Shift (seconds)');
    ylabel('STD positions');
    
    subplot(1,3,3)
    hold on
    title('After correction');
    plot(timestamps,positions);
    % Plot events and corresponding pos. w/ optimal delay
    plti = evidcs - delay*fs;
    plti = plti(plti>0 & plti<=length(positions));
    scatter(timestamps(plti), positions(plti));
end

end