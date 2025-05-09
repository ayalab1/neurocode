function wait_until(targetTime, verbose)
% wait_until Pauses execution until the specified target time.
%
%   wait_until(targetTime) pauses the program until the specified targetTime.
%
%   wait_until(targetTime, verbose) controls whether to display status messages.
%       Set verbose = true to display messages (default: false).
%
%   Example:
%       targetTime = datetime('today') + hours(21); % 9 PM today
%       wait_until(targetTime, true); % Wait with verbose output
%
%   Inputs:
%       targetTime - A datetime object specifying when to resume execution.
%       verbose   - (Optional) Logical flag to enable status messages.

% Set default for verbose if not provided
if nargin < 2
    verbose = false;
end

% Get current time
currentTime = datetime('now', 'TimeZone', targetTime.TimeZone);

% Calculate wait duration in seconds
waitDuration = seconds(targetTime-currentTime);

% Check if the target time is in the future
if waitDuration <= 0
    if verbose
        disp('Target time has already passed. No waiting required.');
    end
    return;
end

% Display wait info if verbose
if verbose
    disp(['Current time: ', char(currentTime)]);
    disp(['Target time: ', char(targetTime)]);
    disp(['Waiting for ', num2str(waitDuration), ' seconds...']);
end

% Pause execution
pause(waitDuration);

% Confirm completion if verbose
if verbose
    disp(['Resuming at: ', char(datetime('now', 'TimeZone', targetTime.TimeZone))]);
end
end