% make an arduino variable to control the signal to the intan
useArduino = false;
if useArduino, aHandle = arduino(); end
pin = 46; % used to be 12 in old code
write = false;

nargin = 0;
% History:
% 3/1/9 mk  Written.

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;
Screen('Preference', 'VisualDebugLevel', 1); % don't show start screen\
Screen('Preference', 'SkipSyncTests', 1);

% create sound ramps:
sound = SoundRamp(20000,5000,0.8,200000,1); soundUpward = sound;
handle_upwardRamp = audioplayer(sound(:,2),1/mode(diff(sound(:,1))));
sound = SoundRamp(5000,20000,0.8,200000,1); soundDownward = sound;
handle_downwardRamp = audioplayer(sound(:,2),1/mode(diff(sound(:,1))));

% Initial stimulus parameters for the grating patch:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = fullfile(pwd,['task-' datestr(clock,30) '-audiovisual-pairing.times']);
stimuli = [0 90]; %stimulus angles: 8 stimuli; -1 is black
nStimuli = length(stimuli);
nRepetitions = 30;
blockSize = 3;
flashDuration = 1;
nFlashesPerStimulus = 2;
nProbeTrials = 5; % for each stimulus, this number of trials will be probe trials (presented sound, but no visual stimulus: this is to check if we get V1 representation of the missing visual cue associated with the presented sound!
firstPossibleProbeTrial = 10; % probe trials will be randomly distributed after the first e.g. 10 trials (it's good if training doesn't start with a probe trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5 || isempty(internalRotation)
    internalRotation = 45;
end

if internalRotation
    rotateMode = kPsychUseTextureMatrixForRotation;
else
    rotateMode = [];
end
if nargin < 4 || isempty(gratingsize)
    gratingsize = 800;
end

% res is the total size of the patch in x- and y- direction, i.e., the
% width and height of the mathematical support:
res = [gratingsize gratingsize];

if nargin < 3 || isempty(freq)
    % Frequency of the grating in cycles per pixel: Here 0.01 cycles per pixel:
    freq = 4.5/gratingsize; %0.05cycle/degree
end

if nargin < 2 || isempty(cyclespersecond)
    cyclespersecond = 1;
end

%if nargin < 1 || isempty(stimuli)
% Tilt stimuli of the grating:
%stimuli = 30;
%end

% Amplitude of the grating in units of absolute display intensity range: A
% setting of 0.5 means that the grating will extend over a range from -0.5
% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
% displayable range. As we select a background color and offset for the
% grating of 0.5 ( == 50% nominal intensity == a nice neutral gray), this
% will extend the sinewaves values from 0 = total black in the minima of
% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
% than 0.5 don't make sense, as parts of the grating would lie outside the
% displayable range for your computers displays:
amplitude = 1;

% Choose screen with maximum id - the secondary display on a dual-display
% setup for display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display, choose a background
% color of 128 = gray, i.e. 50% max intensity:
win = Screen('OpenWindow', screenid, 0);

% Make sure the GLSL shading language is supported:
AssertGLSL;

% Retrieve video redraw interval for later control of our animation timing:
ifi = Screen('GetFlipInterval', win);

% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
phase = 0;

% Compute increment of phase shift per redraw:
phaseincrement = cyclespersecond * 360*0.5;

% Build a procedural sine grating texture for a grating with a support of
% res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
gratingtex = CreateProceduralSineGrating(win, res(1), res(2), [0.5 0.5 0.5 0.0]);

% Wait for release of all keys on keyboard, then sync us to retrace:
KbReleaseWait;

data = nan(1+nRepetitions*nStimuli*(nFlashesPerStimulus+1),3);
% pre-set the pseudorandom stimulus order:
trialType = [];
rng(floor(datenum(clock))); % set a new random number generator every day
for j = 1:nRepetitions/blockSize
    key = 2*ones(blockSize*nStimuli,1); key(1:end/2) = 1;
    indices = randperm(numel(stimuli)*blockSize);
    theseAngles = stimuli(key(indices))';
    trialType = [trialType; theseAngles(:)];
end

% Select random trials to be probe trials
probeTrial = false(size(trialType));
if nProbeTrials>0
    allPossibleLines = 1:nRepetitions*nStimuli;
    % write "-1" in 3rd column for probe trials:
    ok = false;
    while ~ok % start over if you got stuck in an infinite loop
        probeTrial = false(size(trialType));
        for i = 1:length(stimuli)
            if ~(stimuli(i)>-1) continue; end % there is no sound before the control stimulus anyway (or the blank screen (NaN))
            possibleLines = allPossibleLines(trialType==stimuli(i));
            possibleLines(1:firstPossibleProbeTrial-1) = []; % first X trials are off limits and should not be probe trials
            ok = false;
            timestamp = GetSecs;
            while ~ok && GetSecs - timestamp<2 % if you haven't found a solution in 2s, then we're probably stuck and should start over
                indices = sort(randperm(length(possibleLines),nProbeTrials));
                lines = possibleLines(indices)';
                ok = true;
                if any(diff(indices) == 1), ok = false; end % avoid having 2 probe trials in a row (for the same stimulus)
                if any(probeTrial(ismember(allPossibleLines,[lines-1; lines+1]))), ok = false; end % avoid having 2 probe trials in a row (even for different stimuli)
            end
            probeTrial(lines) = -1; % selected for a probe trial
        end
    end
end

startTime = Screen('Flip', win);
data(1,2:3) = [-2 -1]; % start time is code "-2" (for session start) and flash status "-1" (for blank screen)
data(1,1) = datenum(datetime('now','Format','MM/dd/yy HH:mm:ss.SSSSSS')); % save the session start time in the complete format, including the day. This is time 0 and the rest of the timestamps wout be in seconds relative to this start time
currentLine = 2; currentStimulus = 0;
closeScript = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:nRepetitions
    if closeScript
        break;
    end
    for i = 1:length(stimuli)
        if closeScript
            break;
        end
        Screen(win,'FillRect',[0 0 0],[0 0 800 800]);
        timestamp = Screen('Flip', win);
        %     [datestr(clock) ':' num2str(i)]
        if useArduino, for jj = 1:length(pin), writeDigitalPin(aHandle, ['D' num2str(pin(jj))], 0); end; end % send arduino signal 0 (empty screen)

        data(currentLine,1) = timestamp- startTime;
        data(currentLine,2) = -2; % blank screen is coded as stimulus "-2"
        data(currentLine,3) = -1; % flash status "-1" is another code for blank screen

        currentLine = currentLine+1;
        currentStimulus = currentStimulus+1;
        flashStatus = 0;
        randomNumber = rand(1)*10 + 5;
        randomNumber = randomNumber-0.8; randomNumber = 0.2; % for debugging only
        pause(randomNumber);
       
        for k = 1:nFlashesPerStimulus

            if KbCheck
                [keyIsDown, seconds, keyCode ] = KbCheck; % if any keyboard key is pressed
                if keyCode(123) == 0 % if the key was not F12
                    continue; % ignore key press
                end
                closeScript = true; % if F12 was pressed, abort the script
                break
            end

            % Update some grating animation parameters:
            stimulus = trialType(currentStimulus); % pre-set in lines above
            if stimulus<0, actualAmplitude = 0; else, actualAmplitude = amplitude; end % stimulus "-1" is the control stimulus (grey screen, amplitude 0)
            % Increment phase by 180 degrees (to flip it)
            phase = phase + phaseincrement;

            if k == 1
                if stimulus == 0
                    InitializePsychSound(1); nrchannels = 1;
                    pahandle = PsychPortAudio('Open', [], [], 0, round(1/mode(diff(sound(:,1)))), 1);
                    PsychPortAudio('FillBuffer', pahandle, soundUpward(:,2)');
                    timestamp = PsychPortAudio('Start', pahandle, 1, 1, 1);
                elseif stimulus == 90
                    InitializePsychSound(1); nrchannels = 1;
                    pahandle = PsychPortAudio('Open', [], [], 0, round(1/mode(diff(sound(:,1)))), 1);
                    PsychPortAudio('FillBuffer', pahandle, soundDownward(:,2)');
                    timestamp = PsychPortAudio('Start', pahandle, 1, 1, 1);
                end
                data(currentLine,1) = timestamp - startTime;
                data(currentLine,2) = stimulus;
                data(currentLine,3) = 3;
                currentLine = currentLine+1;
                pause(0.8);
            end

            % Draw the grating, centered on the screen, with given r otation 'stimuli',
            % sine grating 'phase' shift and amplitude, rotating via set
            % 'rotateMode'. Note that we pad the last argument with a 4th
            % component, which is 0. This is required, as this argument must be a
            % vector with a number of components that is an integral multiple of 4,
            % i.e. in our case it must have 4 components:
            if ~probeTrial(currentStimulus) % if this is a not a probe trial
                Screen('DrawTexture', win, gratingtex, [], [], stimulus, [], [], [], [], rotateMode, [phase, freq, actualAmplitude, 0]); % draw an empty structure instead of the actual stimulus
            end
            % Show it at next retrace:

            if k == 1
                timestamp = Screen('Flip', win);
            else
                timestamp = Screen('Flip', win);
            end

            if useArduino, for jj = 1:length(pin), writeDigitalPin(aHandle, ['D' num2str(pin(jj))], 1); end; end % send arduino signal
            data(currentLine,1) = timestamp - startTime;
            data(currentLine,2) = stimulus;
            data(currentLine,3) = flashStatus; % flipped compared to the previous one
            if probeTrial(currentStimulus), data(currentLine,3) = -1; end % 2 for blank screen
            currentLine = currentLine+1;
            if write, dlmwrite(filename,data,'precision','%.8f'); end
            flashStatus = 1-flashStatus;
            % Plot progress
            figure(1); clf; plot(trialType,'.-'); ylim(ylim+[-1 1]*30); hold on
            plot(find(probeTrial),trialType(probeTrial),'ro','linewidth',2); ylim(ylim+[-1 1]*30); 
            PlotHVLines(currentStimulus+0.5,'v','k','linewidth',2);
            xlabel('trial ID'); ylabel('orientation (degrees)'); set(gca,'ytick',unique(stimuli)); 
            legend('trial orientations','probe trials','current time');
            legend('box','off');
            set(gca,'box','off','fontsize',12,'TickDir','out');
            drawnow
            pause(flashDuration);
        end
    end
end

Screen(win,'FillRect',[0 0 0],[0 0 800 800]);
timestamp = Screen('Flip', win, timestamp + flashDuration);
pause(2);

% We're done. Close the window. This will also release all other ressources:
sca;
if useArduino, fclose(aHandle); end
beep

% Bye bye!
return;
