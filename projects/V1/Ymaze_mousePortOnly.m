% make an arduino variable to control the signal to the intan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%​
%% Set up parameters​
write = false; % Do we want matlab to save a text file with the timestamps it presented the stimuli? This should be "true" in experiments, but it's useful to disable it (false) while debugging the code
sounds = false; % Do we want to present sounds before each visual cue (inference task)? This option is not yet (today on Jan 28 2022) implemented. !TODO
interruptable = 122; % interruptable: pressing this key should abort the script (key 123 = F12, 122 = F11)​
% The trials are distributed psuedorandomly, in blocks of e.g. 10. Then every 10 trials, we have exactly X left and Y left trials.
blockSize = 10; % This ensures that left/right ratio is maintained throughout the session (rather than all the left trials being assigned randomly to the start of the session, which would introduce a bias in the analyses)
nLeftPerBlock = 5; % the remainder will be right trials (so in every block, we have e.g. 5 left and 5 right trials). Set this to 10 or 0 for training with a single rewarded side.
maxAllowedNumberOfBlocks = 1000; % Trial order is pre-computed at the start of the script. How many blocks do we need to pre-compute to be sure we have enough for the whole task? I think 1000 blocks of 10 trials should be more than enough.
stimuli = [0;90]; % orientation of the two stimuli in degrees (don't change please)
if write, filename = fullfile(pwd,['task-ymaze-' datestr(clock,30) '.times']); end % where to write the matral text file. This is not yet (today on Jan 28 2022) implemented. !TODO
flashDuration = 1; % The stimulus will flip every (e.g. 1) seconds while the animal is blocking the mouse port beam until he moves away
% Set up the minumum time that the animal must break the beam in the mouse port to initiate a new trial:
mimumunHold = 1; % Any shorter time that would not be enough for the animal to see the stimulus. This is not yet (today on Jan 28 2022) implemented. !TODO
maximumHold = 2; % The animal is not allows to get reward at the mouse port forever
arduinoCheckWaitTime = 0.0001;
delay = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Execute preliminary computations

% Arduino stuff:
aHandle = serial('COM6','BaudRate',9600); % load arduino
fopen(aHandle); % open arduino port

% Trial order stuff:
% precompute the blocks (the order of presentation of vertical vs horizontal stimuli​​
[~,blockMatrix] = sort(rand(maxAllowedNumberOfBlocks,blockSize),2); % assign random order (the order of random numbers)
key = 2*ones(1:blockSize); key(1:nLeftPerBlock) = 1; % a "key" vector that translates that number order (from 1 to 10) into "left" and "right" (1=left, 2=right).
blockMatrix = key(blockMatrix); % translate the block matrix using the "key" variable into these [1,2] format
presentationOrder = reshape(blockMatrix',[],1); % the final order is a single vector (we don't care about blocks any more), just "left" vs "right" trials one after the other​

% Sound stuff
if sounds % create sound ramps to play later:
    sound = SoundRamp(20000,5000,0.8,200000,1); soundUpward = sound;
    handle_upwardRamp = audioplayer(sound(:,2),1/mode(diff(sound(:,1))));
    sound = SoundRamp(5000,20000,0.8,200000,1); soundDownward = sound;
    handle_downwardRamp = audioplayer(sound(:,2),1/mode(diff(sound(:,1))));
end

% Psychtoolbox: graphics stuff
AssertOpenGL; % Configure OpenGL Psychtoolbox:
Screen('Preference', 'VisualDebugLevel', 1); % don't show start screen\
Screen('Preference', 'SkipSyncTests', 1); % skip some tests that may trigger code errors​
gratingsize = 800; % in pixels
res = [gratingsize gratingsize]; % res is the total size of the patch in x- and y- direction, i.e., the width and height of the mathematical support:
% Frequency of the grating in cycles per pixel. We want it to match 0.05cycles/degree as in the literature.
freq = 0.05/(4*gratingsize/360); % Since 4 screens cover the whole visual field, then we have 4*800 pixels = 360 degrees, so 1 degree is  4*800/360 pixels
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
screenid = max(Screen('Screens')); screenid = 2; % Choose screen with maximum id - secondary display on a dual-display: change this to 1,2 or 3 if the stimuli appear on the wrong screen
win = Screen('OpenWindow', screenid, 0); % Open a fullscreen onscreen window on that display
AssertGLSL; % Make sure the GLSL shading language is supported:​
gratingtex = CreateProceduralSineGrating(win, res(1), res(2), [0.5 0.5 0.5 0.0]); % Build a grating texture with res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
KbReleaseWait; % Wait for release of all keys on keyboard, then sync us to retrace (not sure if that's needed). Test without it !TODO​
timestamp = Screen('Flip', win); % collect Psychtoolbox's start time for the text file records !TODO
data(1,1) = datenum(datetime('now','Format','MM/dd/yy HH:mm:ss.SSSSSS'));
data(1,[2 3]) = [-1 2]; % record start time with label (-1 2) to differentiate it from other timestamps should be a bit different from the rest
currentLine = 2; % for writing the data matrix later. !TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the loop for the task
currentStimulus = 0; % This variable tells us which row of "presentationOrder" are we currently on (we haven't started yet so it's 0)
mode = 1; % mode 1 = mouse port, mode 2 = left, mode 3 = reward, mode 0 = timeout, mode 4=readdata, mode 5=report screen change. Start with mode 1.
fprintf(aHandle,'%s',char(mode)); % write command to arduino
abortScript = false; % This variable keeps track of when we should abort the experiment (aborted e.g. by pressing F12)
variable = []; vari = []; tic;
flushinput(aHandle); flushoutput(aHandle);
while ~abortScript
    display(['In mode ' num2str(mode)])
    switch mode % check the mode and execute the corresponding code
        case 1 % mouse port
            if abortScript % make sure we exit if experimenter aborted the experiment
                break
            end

            Screen(win,'FillRect',[0 0 0],[0 0 800 800]); % Start with a blank screen:
            timestamp = Screen('Flip', win); % write code to save timestamp later !TODO

            mousePortCrossed = false;
            

            % Prepare the stimulus for when the mouse port beam gets crossed
            currentStimulus = currentStimulus+1;
            direction = presentationOrder(currentStimulus);
            stimulus = stimuli(direction);
            phase = 0; % Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
            Screen('DrawTexture', win, gratingtex, [], [], stimulus, [], [], [], [], rotateMode, [phase, freq, amplitude, 0]); % the stimulus is prepared but will not appear on screen until we send the Flip command

            % Wait for the mouse to approach the mouse port and block the beam:
            while ~mousePortCrossed
                % In this loop, we will check with the arduino if the mouse has crossed the mouse port beam
                % and also check if the experimenter pressed F12 to abort the experiment
                vari(end+1,3) = toc;
                fwrite(aHandle,4,'uchar');% ask arduino send info
                vari(end,[1 2]) = [aHandle.BytesAvailable toc];
                if(vari(end,1)>6) % receive it if there is one
                    signal0 =  fread(aHandle,aHandle.BytesAvailable);
                    signalEnd = signal0(end-6:end); % the last 7 digits contain 6 with the signal and 1 with the marker "250"
                    signal = circshift(signalEnd,-find(signalEnd==250)); % shift them so that the marker "250" is last
                    mousePortCrossed = signal(1);
                    leftPortCrossed = signal(2);
                    rightPortCrossed = signal(3);
                    variable(end+1,:) = signal;
                    for i=1:6, subplot(3,2,i); plot(variable(:,i)); title(variable(end,i)); end; drawnow
                    flushinput(aHandle); flushoutput(aHandle);
                end
                if KbCheck
                    [keyIsDown, seconds, keyCode ] = KbCheck; % if any keyboard key is pressed
                    if keyCode(interruptable)==0 % if the key was not F12
                        continue; % ignore key press
                    end
                    abortScript = true; % if F12 was pressed, abort the script
                    break
                end
                pause(arduinoCheckWaitTime); % check the beam or if experimenter aborted the experiment every 10 ms
            end

            timestamp = GetSecs; % Record the timestamp arduino reported that the animal crossed! This will help with synchronization.
            if abortScript % make sure we exit if experimenter aborted the experiment
                break
            end
            % If we exited the while loop without aborting the script, this means that the mouse port beam is currently being crossed
%             timestamp = Screen('Flip', win); % Make the stimulus appear! Raly: write code to save timestamp later
            timestamp0 = timestamp;
            while mousePortCrossed && ~abortScript && GetSecs - timestamp0 < maximumHold
                vari(end+1,3) = toc;
                fwrite(aHandle,4,'uchar');% ask arduino send info
                vari(end,[1 2]) = [aHandle.BytesAvailable toc];
                if(aHandle.BytesAvailable>6) % receive it if there is one
                    signal0 =  fread(aHandle,aHandle.BytesAvailable);
                    signalEnd = signal0(end-6:end); % the last 7 digits contain 6 with the signal and 1 with the marker "250"
                    signal = circshift(signalEnd,-find(signalEnd==250)); % shift them so that the marker "250" is last
                    mousePortCrossed = signal(1);
                    leftPortCrossed = signal(2);
                    rightPortCrossed = signal(3);
                    variable(end+1,:) = signal;
                    for i=1:6, subplot(3,2,i); plot(variable(:,i)); title(variable(end,i)); end; drawnow
                    flushinput(aHandle); flushoutput(aHandle);
                end
                if KbCheck
                    [keyIsDown, seconds, keyCode ] = KbCheck; % if any keyboard key is pressed
                    if keyCode(123)==0 % if the key was not F12
                        continue; % ignore key press
                    end
                    abortScript = true; % if F12 was pressed, abort the script
                    fprintf(aHandle,'%s',char(0)); % write command to arduino to stop all reward ports
                    break
                end

                if GetSecs - timestamp > 1 % if more than 1s has passed since the stimulus has been shown (or flipped), flip the screen
                    phase = phase + 180;
                    Screen('DrawTexture', win, gratingtex, [], [], stimulus, [], [], [], [], rotateMode, [phase, freq, amplitude, 0]);
%                     timestamp = Screen('Flip', win); % Raly: write code to save timestamp later
                end
                pause(arduinoCheckWaitTime); % check the beam or if experimenter aborted the experiment every 10 ms

            end
            % When the mouse moves away from the mouse port, put a blank screen immediately
            % (the mouse can view stimuli only while it's blocking the mouse port beam)
            Screen(win,'FillRect',[0 0 0],[0 0 800 800]);
            timestamp = Screen('Flip', win); % Raly: write code to save timestamp later !TODO
            fprintf(aHandle,'%s',char(5)); % tell arduino to report screen change to intan
            mode = 1; % change the mode. Reward should now be available at the Left or the Right reward port
            pause(delay)
            fprintf(aHandle,'%s',char(mode)); % change mode: reward should now be at the left/right arm​
            flushinput(aHandle); flushoutput(aHandle);
%         case 2
%             pause(delay)
%             mode = 1;
    end
end

% We're done. Close the window. This will also release all other ressources:
sca;
fprintf(aHandle,'%s',char(0)); % write command to arduino to stop all reward ports
fclose(aHandle);
return