function wait_until(targetTime, verbose, showCountdown)
% wait_until Pauses execution until the specified target time.
%
%   wait_until(targetTime) pauses the program until the specified targetTime.
%
%   wait_until(targetTime, verbose) controls whether to display status messages.
%       Set verbose = true to display messages (default: false).
%
%   wait_until(targetTime, verbose, showCountdown) enables a live countdown.
%       Set showCountdown = true to display a dynamic countdown (default: false).
%
%   Example:
%       targetTime = datetime('today') + hours(21); % 9 PM today
%       wait_until(targetTime, true, true); % Wait with verbose output + countdown
%
%   Inputs:
%       targetTime    - A datetime object specifying when to resume execution.
%       verbose      - (Optional) Logical flag to enable status messages.
%       showCountdown - (Optional) Logical flag to enable live countdown.

% Set defaults
if nargin < 2
    verbose = false;
end
if nargin < 3
    showCountdown = false;
end

% Check for active breakpoints
dbstatus_ = dbstatus('-completenames');
if ~isempty(dbstatus_)
    warning(['Active breakpoints detected! ', ...
        'Consider removing breakpoints (dbclear all) before using wait_until.']);
end

% Get current time (match timezone of targetTime)
currentTime = datetime('now', 'TimeZone', targetTime.TimeZone);
waitDuration = seconds(targetTime-currentTime);

% Check if target time is in the past
if waitDuration <= 0
    if verbose
        disp('Target time has already passed. No waiting required.');
    end
    return;
end

% Display initial info if verbose
if verbose
    disp(['Current time: ', char(currentTime)]);
    disp(['Target time: ', char(targetTime)]);
    disp(['Waiting for ', num2str(waitDuration), ' seconds...']);
end

% If countdown is enabled, display updates every second
if showCountdown
    fprintf('Count: ');
    remainingTime = waitDuration;
    while remainingTime > 0
        % Format remaining time as HH:MM:SS
        hoursLeft = floor(remainingTime/3600);
        minutesLeft = floor(mod(remainingTime, 3600)/60);
        secondsLeft = floor(mod(remainingTime, 60));
        countdownStr = sprintf('%02d:%02d:%02d', hoursLeft, minutesLeft, secondsLeft);

        % Overwrite previous countdown line
        fprintf([repmat('\b', 1, 8), '%s'], countdownStr);

        % Wait 1 second
        pause(1);
        remainingTime = remainingTime - 1;
    end
    fprintf('\n'); % Move to next line after countdown
else
    % If no countdown, just pause for the full duration
    pause(waitDuration);
end

% Display completion message if verbose
if verbose
    disp(['Resuming at: ', char(datetime('now', 'TimeZone', targetTime.TimeZone))]);
end
end