function [data,arduinoOutput] = Ymaze(varargin)

%Ymaze - Run the Ymaze in AyAlab!
%
% Select the parameters you want and start the recording.
%
%  USAGE
%
%    [data,arduinoOutput] = Ymaze(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'write'           boolean indicating if timestamps should be saved on disk (default = true)
%                       it may be useful to disable it (set it to "false") while debugging the code
%     'folder'          folder where results will be saved (default = current folder)
%     'interruptable'   keyboard key which will interrupt the script (default = 123, i.e. key F12)
%     'screenid'        the id of the screen where stimuli should be shown (default = 2)
%     'blockSize'       number of trials composing a block (default = 4). The trials are distributed pseudorandomly: 
%                       by default, every block of 10 trials, we have exactly 2 left and 2 trials. 
%                       This ensures that left/right ratio is maintained throughout the session.
%     'nLeftPerBlock'   number of "left" trials in a block (default = 2). When set to 0, all blocks will be
%                       composed exclusively of right trials (and the reverse for nLeftPerBlock=blockSize).
%     'flashDuration'   time before stimulus flips (default = 1 second) as the animal is blocking the mouse port beam
%     'minimumHold'     minumum time that the animal must break the beam in the mouse port to initiate a new trial
%                       (default = 1 second). Shorter times (not enough for the animal to see the stimulus) would 
%                       result in repeat trials (animal must nose poke again).
%     'maximumHold'     maximim time that the animal would get rewards at the mouse port (default = 2 seconds).
%                       Past this limit, a new trial is initiated and reward is no longer available at the mouse port.
%     'waitBetweenBeamAndReward'     time waited before dispensing water from the left/right reward ports (after the animal crossed the beam).
%                       This is to allow the animal to drink the reward as it is dispensed.
%     'show'            boolean indicating if arduino input and task progress should be displayed (default = true)
%     'sounds'          boolean indicating if sounds should precede cue presenation (default = false)
%                       NOTE: this is not yet implemented!!! TODO
%    =========================================================================
%
%  OUTPUT
%
%    data               final matrix of timestamps and events that matlab will save for future synchronization
%    arduinoOutput      exact readings of the arduino variables (beam crossings and number of rewards)
%
%
% Copyright (C) 2022 by Ralitsa Todorova and Hongyu Chang
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%​
%% Set up default parameters

folder = pwd;
write = true;
interruptable = 123; % F12 key
screenid = 2;
blockSize = 4;
nLeftPerBlock = 2;
flashDuration = 1;
minimumHold = 0.5;
maximumHold = 5;
waitBetweenBeamAndReward = 1;
show = true;
sounds = false; % Not yet implemented! TODO

% Parameters that are not set by the user:
% Trial order is pre-computed at the start of the script. How many blocks do we need to pre-compute to be sure we have enough for the whole task? I think 1000 blocks of 10 trials should be more than enough.
maxAllowedNumberOfBlocks = 1000;
stimuli = [0 90]; % orientation of the two stimuli in degrees (don't change please)

% Parse parameter list
for i = 1:2:length(varargin)
	if ~ischar(varargin{i})
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).']);
	end
    switch(lower(varargin{i}))
        case 'write'
            write = varargin{i+1};
            if ~islogical(write)
                error('Incorrect value for property ''write'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'folder'
            folder = varargin{i+1};
            if ~isfolder(folder)
                error('Incorrect value for property ''folder'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'interruptable'
            interruptable = varargin{i+1};
            if ~isnumeric(interruptable)
                error('Incorrect value for property ''interruptable'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'screenid'
            screenid = varargin{i+1};
            if ~isnumeric(screenid)
                error('Incorrect value for property ''screen'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'blocksize'
            blockSize = varargin{i+1};
            if ~isnumeric(blockSize)
                error('Incorrect value for property ''blockSize'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'nleftperblock'
            nLeftPerBlock = varargin{i+1};
            if ~isnumeric(nLeftPerBlock)
                error('Incorrect value for property ''nLeftPerBlock'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'flashduration'
            flashDuration = varargin{i+1};
            if ~isnumeric(flashDuration)
                error('Incorrect value for property ''flashDuration'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'minimumhold'
            minimumHold = varargin{i+1};
            if ~isnumeric(minimumHold)
                error('Incorrect value for property ''minimumHold'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'maximumhold'
            maximumHold = varargin{i+1};
            if ~isnumeric(maximumHold)
                error('Incorrect value for property ''maximumHold'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'waitbetweenbeamrndreward'
            waitBetweenBeamAndReward = varargin{i+1};
            if ~isnumeric(waitBetweenBeamAndReward)
                error('Incorrect value for property ''waitBetweenBeamAndReward'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end
        case 'show'
            show = varargin{i+1};
            if ischar(show), if strcmp(show,'on'), show = true; elseif strcmp(show,'off'), show = false; end; end
            if ~islogical(show)
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end	
     case 'sounds'
            sounds = varargin{i+1};
            if ~islogical(sounds)
                error('Incorrect value for property ''sounds'' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).');
            end

        otherwise
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Ymaze">Ymaze</a>'' for details).']);
	end
end

if nLeftPerBlock== Inf, nLeftPerBlock = blockSize; end
if write, filename = fullfile(folder, ['task-ymaze-' datestr(clock, 30) '.times']); end % where to write the matral text file. This is not yet (today on Jan 28 2022) implemented. !TODO
if write, filename2 = fullfile(folder, ['arduino-outputs-' filename(end-20:end)]); end % where to write the matral text file. This is not yet (today on Jan 28 2022) implemented. !TODO
arduinoCheckWaitTime = 0.0001;

display(['Blocks will be ' num2str(blockSize) ' trials, ' num2str(nLeftPerBlock) ' left and ' num2str(blockSize-nLeftPerBlock) ' right trials']);
if write
    [~,name1] = fileparts(filename); [~,name2] = fileparts(filename2);
    display(['Matlab files will be saved in ' folder ', ' name1 '.times and ' name2 '.times.']);
else
    display(['Matlab files not be saved.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Execute preliminary computations

% Arduino stuff:
aHandle = serial('COM6', 'BaudRate', 9600); % load arduino
fopen(aHandle); % open arduino port

% Trial order stuff:
% precompute the blocks (the order of presentation of vertical vs horizontal stimuli​​
[~, blockMatrix] = sort(rand(maxAllowedNumberOfBlocks, blockSize), 2); % assign random order (the order of random numbers)
key = 2*ones(1:blockSize); key(1:nLeftPerBlock) = 1; % a "key" vector that translates that number order (from 1 to 10) into "left" and "right" (1=left, 2=right).
blockMatrix = key(blockMatrix); % translate the block matrix using the "key" arduinoOutputs into these [1, 2] format
presentationOrder = reshape(blockMatrix', [], 1); % the final order is a single vector (we don't care about blocks any more), just "left" vs "right" trials one after the other​

% Sound stuff
if sounds % create sound ramps to play later:
    sound = SoundRamp(20000, 5000, 0.8, 200000, 1); soundUpward = sound;
    handle_upwardRamp = audioplayer(sound(:, 2), 1/mode(diff(sound(:, 1))));
    sound = SoundRamp(5000, 20000, 0.8, 200000, 1); soundDownward = sound;
    handle_downwardRamp = audioplayer(sound(:, 2), 1/mode(diff(sound(:, 1))));
end

% Psychtoolbox: graphics stuff
AssertOpenGL; % Configure OpenGL Psychtoolbox:
Screen('Preference', 'Verbosity', 0)
Screen('Preference', 'VisualDebugLevel', 0); % don't show start screen\
Screen('Preference', 'SkipSyncTests', 1); % skip some tests that may trigger code errors​
gratingsize = 800; % in pixels
res = [gratingsize gratingsize]; % res is the total size of the patch in x- and y- direction, i.e., the width and height of the mathematical support:
% Frequency of the grating in cycles per pixel. We want it to match 0.05cycles/degree as in the literature.
freq = 0.05/(4*gratingsize/360); % Since 4 screens cover the whole visual field, then we have 4*800 pixels = 360 degrees, so 1 degree is 4*800/360 pixels
% Amplitude of the grating in units of absolute display intensity range: A
% setting of 0.5 means that the grating will extend over a range from -0.5
% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
% displayable range. As we select a background color and offset for the
% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
% will extend the sinewaves values from 0 = total black in the minima of
% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
% than 0.5 don't make sense, as parts of the grating would lie outside the
% displayable range for your computers displays:
amplitude = 1; % We use 1 despite explanation above, because this gives us black clear and white stripes in between the gratings
% We were meant to do tests about which amplitude maximizes the V1
% response, but we haven't done those. Maybe we could, as part of the
% pilot? !TODO​
rotateMode = kPsychUseTextureMatrixForRotation;
% This is now a parameter: screenid = max(Screen('Screens')); screenid = 2; % Choose screen with maximum id - secondary display on a dual-display: change this to 1, 2 or 3 if the stimuli appear on the wrong screen
win = Screen('OpenWindow', screenid, 0); % Open a fullscreen onscreen window on that display
AssertGLSL; % Make sure the GLSL shading language is supported:​
gratingtex = CreateProceduralSineGrating(win, res(1), res(2), [0.5 0.5 0.5 0.0]); % Build a grating texture with res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
KbReleaseWait; % Wait for release of all keys on keyboard, then sync us to retrace (not sure if that's needed). Test without it !TODO​

startTime = Screen('Flip', win);
data(1, 2:3) = [-2 -1]; % start time is code "-2" (for session start) and flash status "-1" (for blank screen)
data(1, 1) = datenum(datetime('now', 'Format', 'MM/dd/yy HH:mm:ss.SSSSSS')); % save the session start time in the complete format, including the day. This is time 0 and the rest of the timestamps wout be in seconds relative to this start time
currentLine = 2;

% A (complete?) guide on the reading the "data" matrix:
% The first column is timestamps. The very first timestamp is computer time of the start of the experiment. Subsequent timestamps are [seconds]
% following the start of the experiment. In older versions (pre- Feb 7 2022), all timestamps are in computer time [days]
% The second column is stimulus code. Values from 0 to 360 indicate degrees of the stimulus that has just appeared. -1 indicates a grey screen
% (control, same luminance as shown stimuli but zero amplitude), -2 indicates a the appearance of a blank screen (completely black), while
% -3 indicates that the screen continues to be blank (no change from before) as the line marks a non-screen event. -100 indicates that the experiment
% has been aborted through pressing the interruption key (usually F12).  Values of -10, -20, and -30 indicates beam crossings to the mouse port, the 
% left port, and the right port, respectively.
% The third column indicates flash status. Values of 0 and 1 indicate the actual flash status of a visual cue, with 0 being the initial version of the
% grating and 1 indicating that it has flashed or flipped phase by 180 degrees (black becomes white and white becomes black). A value of 3 indicates
% that a sound was played (the sound corresponding to the cue indicated in the second column). A value of -1 indicates the apparition of a blank
% screen. A value of -50 indicates an error (crossing the wrong left/right beam), while a value of 50 indicates a correct choice (crossing the correct
% left/right beam).
% The fourth column indicates the current trial number.
% The fifth column indicates the current mode of the task, with a value of 1 indicating that the mouse should direct itself towards the mouse port to
% initiate a trial, 2 indicating that the mouse should go to the left arm for a reward, and 2 indicating that the mouse should go to the right arm for
% a reward. 0 indicates the end of the experiment.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the loop for the task
trialNumber = 0; % This arduinoOutputs tells us which row of "presentationOrder" are we currently on (we haven't started yet so it's 0)
mode = 1; % mode 1 = mouse port, mode 2 = left, mode 3 = reward, mode 0 = timeout, mode 4=readdata, mode 5=report screen change. Start with mode 1.
fprintf(aHandle, '%s', char(mode)); % write command to arduino
abortScript = false; % This arduinoOutputs keeps track of when we should abort the experiment (aborted e.g. by pressing F12)
arduinoOutputs = []; arduinoOutputs(1,:) = [data(1) zeros(1,6) 0];
flushinput(aHandle); flushoutput(aHandle);
timeLeft = minimumHold; repeatedTrial = false;
mousePortCrossed = false; leftPortCrossed = false; rightPortCrossed = false;
while ~abortScript
    display(['In mode ' num2str(mode)]);
    switch mode % check the mode and execute the corresponding code
        case 1 % mouse port
            if abortScript, break; end % make sure we exit if experimenter aborted the experiment

            Screen(win, 'FillRect', [0 0 0], [0 0 800 800]); % Start with a blank screen:
            timestamp = Screen('Flip', win); %
            mousePortCrossed = false;
            data(currentLine, 1) = timestamp - startTime;
            data(currentLine, 2) = -3; % no change is coded by "-3"
            data(currentLine, 3) = -1; % flash status "-1" is another code for blank screen
            data(currentLine, 4) = trialNumber; % current trial number
            data(currentLine, 5) = mode; % current mode (1=mouse port, 2=left port, 3=right port, 5=timeout)
            currentLine = currentLine+1;
            if write, dlmwrite(filename, data, 'precision', '%.8f'); dlmwrite(filename2, arduinoOutputs, 'precision', '%.8f'); end

            % Prepare the stimulus for when the mouse port beam gets crossed
            if ~repeatedTrial, trialNumber = trialNumber+1; end % go to next trial, unless we are repeating a trial in which the animal didn't stay at the mouse port the minimum amount of time
            direction = presentationOrder(trialNumber);
            stimulus = stimuli(direction);
            phase = 0; % Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
            Screen('DrawTexture', win, gratingtex, [], [], stimulus, [], [], [], [], rotateMode, [phase, freq, amplitude, 0]); % the stimulus is prepared but will not appear on screen until we send the Flip command

            % Wait for the mouse to approach the mouse port and block the beam:
            while ~mousePortCrossed
                % In this loop, we will check with the arduino if the mouse has crossed the mouse port beam
                % and also check if the experimenter pressed F12 to abort the experiment
                checkArduino; % this script is found in the "private" folder
                checkInterruptionKey; % this script is found in the "private" folder
                if abortScript, break; end
                pause(arduinoCheckWaitTime); % check the beam or if experimenter aborted the experiment every X ms
            end
            if abortScript, break; end % make sure we exit if experimenter aborted the experiment

            % If we exited the while loop without aborting the script, this means that the mouse port beam is currently being crossed
            timestamp = Screen('Flip', win); % Make the stimulus appear!
            data(currentLine, 1) = timestamp - startTime;
            data(currentLine, 2) = stimulus; % right port crossing is coded as stimulus "-30"
            data(currentLine, 3) = 0; % flash status "0" is the initial orientation (not flipped like during a flash)
            data(currentLine, 4) = trialNumber; % current trial number
            data(currentLine, 5) = mode; % current mode (1=mouse port, 2=left port, 3=right port, 5=timeout)
            if write, dlmwrite(filename, data, 'precision', '%.8f'); dlmwrite(filename2, arduinoOutputs, 'precision', '%.8f'); end
            currentLine = currentLine + 1;
            timestamp0 = timestamp; flashTimestamp = timestamp;
            while mousePortCrossed && ~abortScript && (GetSecs - timestamp0 < maximumHold)
                checkArduino
                checkInterruptionKey; % this script is found in the "private" folder
                if abortScript, break; end
                if GetSecs - flashTimestamp > flashDuration % if more than 1s has passed since the stimulus has been shown (or flipped), flip the screen
                    phase = phase + 180;
                    Screen('DrawTexture', win, gratingtex, [], [], stimulus, [], [], [], [], rotateMode, [phase, freq, amplitude, 0]);
                    flashTimestamp = Screen('Flip', win); % Raly: write code to save timestamp later
                    data(currentLine, 1) = timestamp - startTime;
                    data(currentLine, 2) = stimulus; % right port crossing is coded as stimulus "-30"
                    data(currentLine, 3) = phase~=0; % flash status "1" is the the flipped version of the original orientation
                    data(currentLine, 4) = trialNumber; % current trial number
                    data(currentLine, 5) = mode; % current mode (1=mouse port, 2=left port, 3=right port, 5=timeout)
                    currentLine = currentLine + 1;
                end
                pause(arduinoCheckWaitTime); % check the beam or if experimenter aborted the experiment every 10 ms
            end
            timeLeft = 0; % don't need this functionality in this training version where the cue stays on during the whole trial
            timestamp = GetSecs;
            if timeLeft<=0 % if animal didn't stay at the mouse port the minimum time, don't change the mode and repeat the mouse port mode (1)
                fprintf(aHandle, '%s', char(5)); % tell arduino to report screen change to intan
                mode = direction+1; % change the mode. Reward should now be available at the Left or the Right reward port
                fprintf(aHandle, '%s', char(mode)); % change mode: reward should now be at the left/right arm​
                flushinput(aHandle); flushoutput(aHandle);
                timeLeft = minimumHold; % reset the time left for the next stimulus
                repeatedTrial = false;
            else
                repeatedTrial = true;
            end
        case 2 % Left
            if abortScript, break; end % make sure we exit if experimenter aborted the experiment

            % Wait for the mouse to approach the mouse port and block the beam:
            leftPortCrossed = false; rightPortCrossed = false;
            while ~leftPortCrossed && ~rightPortCrossed
                % In this loop, we will check with the arduino if the mouse has crossed any beam
                % and also check if the experimenter pressed F12 to abort the experiment
                % In this loop, we will check with the arduino if the mouse has crossed the mouse port beam
                % and also check if the experimenter pressed F12 to abort the experiment
                checkArduino; % this script is found in the "private" folder
                checkInterruptionKey; % this script is found in the "private" folder
                if abortScript, break; end
                pause(arduinoCheckWaitTime); % check the beam or if experimenter aborted the experiment every X ms
            end

            % If we exited the while loop without aborting the script, this
            % means that one of the reward beams was crossed
            if leftPortCrossed % Correct!
                pause(waitBetweenBeamAndReward) % wait for five seconds for the animal to consume its reward before changing the mode !TODO decide how long the wait should be
                mode = 1; % reward should now be back at the mouse port
                fprintf(aHandle, '%s', char(mode)); % write command to arduino
                flushinput(aHandle); flushoutput(aHandle);
                display('Correct!')
                data(currentLine-1, 3) = 50; % flash status "50" is indicates reward
            elseif rightPortCrossed % Error!
                % No reward. We should think about introducing an error
                % sound or a timeout? For the time being, the animal just
                % needs to go back to the start
                mode = 1; % reward should now be back at the mouse port
                fprintf(aHandle, '%s', char(mode)); % write command to arduino
                flushinput(aHandle); flushoutput(aHandle);
                display('Error!')
                data(currentLine-1, 3) = -50; % flash status "-50" is indicates an error
            end
            % Make the screen blank again
            Screen(win, 'FillRect', [0 0 0], [0 0 800 800]);
            timestamp = Screen('Flip', win);
            data(currentLine, 1) = timestamp - startTime;
            data(currentLine, 2) = -2; % blank screen is coded by "-2"
            data(currentLine, 3) = -1; % flash status "-1" is another code for blank screen
            data(currentLine, 4) = trialNumber; % current trial number
            data(currentLine, 5) = mode; % current mode (1=mouse port, 2=left port, 3=right port, 5=timeout)
            currentLine = currentLine+1;
            if write, dlmwrite(filename, data, 'precision', '%.8f'); end % here we are adding the "correct/error" status to the previously saved arduino signal indicating the beam cross
        case 3 % Right
            if abortScript, break; end % make sure we exit if experimenter aborted the experiment

            % Wait for the mouse to approach the mouse port and block the beam:
            leftPortCrossed = false; rightPortCrossed = false;
            while ~leftPortCrossed && ~rightPortCrossed
                % In this loop, we will check with the arduino if the mouse has crossed the mouse port beam
                % and also check if the experimenter pressed F12 to abort the experiment
                checkArduino; % this script is found in the "private" folder
                checkInterruptionKey; % this script is found in the "private" folder
                if abortScript, break; end
                pause(arduinoCheckWaitTime); % check the beam or if experimenter aborted the experiment every X ms
            end
            % If we exited the while loop without aborting the script, this
            % means that one of the reward beams was crossed
            if rightPortCrossed% Correct!
                pause(waitBetweenBeamAndReward) % wait for five seconds for the animal to consume its reward before changing the mode !TODO decide how long the wait should be
                mode = 1; % reward should now be back at the mouse port
                fprintf(aHandle, '%s', char(mode)); % write command to arduino
                flushinput(aHandle); flushoutput(aHandle);
                display('Correct!')
                data(currentLine-1, 3) = 50; % flash status "50" is indicates reward
            elseif leftPortCrossed % Error!
                % No reward. We should think about introducing an error
                % sound or a timeout? For the time being, the animal just
                % needs to go back to the start
                mode = 1; % reward should now be back at the mouse port
                fprintf(aHandle, '%s', char(mode)); % write command to arduino
                flushinput(aHandle); flushoutput(aHandle);
                display('Error!')
                data(currentLine-1, 3) = -50; % flash status "-50" is indicates an error
            end
            % Make the screen blank again
            Screen(win, 'FillRect', [0 0 0], [0 0 800 800]);
            timestamp = Screen('Flip', win);
            data(currentLine, 1) = timestamp - startTime;
            data(currentLine, 2) = -2; % blank screen is coded by "-2"
            data(currentLine, 3) = -1; % flash status "-1" is another code for blank screen
            data(currentLine, 4) = trialNumber; % current trial number
            data(currentLine, 5) = mode; % current mode (1=mouse port, 2=left port, 3=right port, 5=timeout)
            currentLine = currentLine+1;
            if write, dlmwrite(filename, data, 'precision', '%.8f'); dlmwrite(filename2, arduinoOutputs, 'precision', '%.8f'); end

            if write, dlmwrite(filename, data, 'precision', '%.8f'); end % here we are adding the "correct/error" status to the previously saved arduino signal indicating the beam cross
    end
end

% We're done. Close the window. This will also release all other ressources:
sca;
fprintf(aHandle, '%s', char(0)); % write command to arduino to stop all reward ports
fclose(aHandle);
if write, dlmwrite(filename, data, 'precision', '%.8f'); dlmwrite(filename2, arduinoOutputs, 'precision', '%.8f'); end
return
