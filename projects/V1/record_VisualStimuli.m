% function task(cyclespersecond, freq, gratingsize, internalRotation)
% function DriftDemo4([stimuli=0][, cyclespersecond=1][, freq=1/360][, gratingsize=360][, internalRotation=0])
% ___________________________________________________________________
%
% Display an animated grating, using the new Screen('DrawTexture') command.
% This demo demonstrates fast drawing of such a grating via use of procedural
% texture mapping. It only works on hardware with support for the GLSL
% shading language, vertex- and fragmentshaders. The demo ends if you press
% any key on the keyboard.
%
% The grating is not encoded into a texture, but instead a little algorithm - a
% procedural texture shader - is executed on the graphics processor (GPU)
% to compute the grating on-the-fly during drawing.clear
%
% This is very fast and efficient! All parameters of the grating can be
% changed dynamically. For a similar approach wrt. Gabors, check out
% ProceduralGaborDemo. For an extremely fast aproach for drawing many Gabor
% patches at once, check out ProceduralGarboriumDemo. That demo could be
% easily customized to draw many sine gratings by mixing code from that
% demo with setup code from this demo.
%
% Optional Parameters:
% 'stimuli' = Rotation stimuli of grating in degrees.
% 'internalRotation' = Shall the rectangular image patch be rotated
% (default), or the grating within the rectangular patch?
% gratingsize = Size of 2D grating patch in pixels.
% freq = Frequency of sine grating in cycles per pixel.
% cyclespersecond = Drift speed in cycles per second.
%

% make an arduino variable to control the signal to the intan
useArduino = true;
if useArduino, a = arduino(); end
write = true;
nargin = 0;
% Choose screen with maximum id - e.g. the secondary display on a dual-display setup for display:
screenid = max(Screen('Screens')); screenid=3;
nRepetitions = 30; 

% History:
% 3/1/9  mk   Written.

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;
Screen('Preference', 'VisualDebugLevel', 1); % don't show start screen\
Screen('Preference', 'SkipSyncTests', 1); 
% Initial stimulus parameters for the grating patch:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = fullfile(pwd,['task-' datestr(clock,30) '.times']);
movieDurationSecs = 2; % how long the display shows each stimulus
stimuli = [-1 0 22.5 45  67.5  90  112.5  135 157.5]; %stimulus angles: 8 stimuli; -1 is black
nStimuli = length(stimuli);
flashDuration = 1;
nFlashesPerStimulus = floor(movieDurationSecs/flashDuration);
flashDuration = movieDurationSecs/nFlashesPerStimulus;
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
% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
% will extend the sinewaves values from 0 = total black in the minima of
% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
% than 0.5 don't make sense, as parts of the grating would lie outside the
% displayable range for your computers displays:
amplitude = 1;

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
phaseincrement =cyclespersecond * 360*0.5;

% Build a procedural sine grating texture for a grating with a support of
% res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
gratingtex = CreateProceduralSineGrating(win, res(1), res(2), [0.5 0.5 0.5 0.0]);

% Wait for release of all keys on keyboard, then sync us to retrace:
KbReleaseWait;
data = nan(1+nRepetitions*nStimuli*(nFlashesPerStimulus+1),3);
% pre-set the pseudorandom stimulus order:
angles = [];
for j=1:nRepetitions
    theseAngles = stimuli(randperm(numel(stimuli)));
    angles = [angles; theseAngles(:)];
end
angles = repelem(angles,nFlashesPerStimulus+1); % add one extra line for the blank screen
ok = diff([-1;angles])~=0; angles(ok) = nan; % first line for each stimulus is actually the line for the blank screen
data(2:end,2) = angles;
vbl = Screen('Flip', win);
data(1,1) = datenum(datetime('now','Format','MM/dd/yy HH:mm:ss.SSSSSS'));
data(1,[2 3]) = [-1 2]; % start time should be a bit different from the rest
vblendtime = vbl + movieDurationSecs;
currentLine = 2;
closeScript = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:nRepetitions
    if closeScript
        break;
    end
    stimuli(:) = stimuli(randperm(numel(stimuli)));
    for i=1:length(stimuli)
        if closeScript
            break;
        end
        
        randomNumber = rand(1)*2+1;
        Screen(win,'FillRect',[0 0 0],[0 0 800 800]);
        vbl = Screen('Flip', win, vbl+flashDuration);
        %         [datestr(clock) ':' num2str(i)]
        if useArduino, writeDigitalPin(a, 'D46', 0);  end % send arduino signal 0 (empty screen)
        data(currentLine,1) = datenum(datetime('now','Format','MM/dd/yy HH:mm:ss.SSSSSS'));
        data(currentLine,2)= -2;
        data(currentLine,3) = 2; % flipped compared to the previous one
        
        currentLine = currentLine+1;
        flashStatus = 0;
        currentStimulusStartTime = vbl + randomNumber;
        for k=1:nFlashesPerStimulus
            if KbCheck
                    [keyIsDown, seconds, keyCode ] = KbCheck; % if any keyboard key is pressed
                    if keyCode(123)==0 % if the key was not F12
                        continue; % ignore key press
                    end
                    closeScript = true; % if F12 was pressed, abort the script
                    break
            end
            % Update some grating animation parameters:
            stimulus = data(currentLine,2); % pre-set in lines above
            if stimulus<0, actualAmplitude = 0; else, actualAmplitude = amplitude; end
            % Increment phase by 1 degree:
            phase = phase + phaseincrement;
            
            % Draw the grating, centered on the screen, with given r otation 'stimuli',
            % sine grating 'phase' shift and amplitude, rotating via set
            % 'rotateMode'. Note that we pad the last argument with a 4th
            % component, which is 0. This is required, as this argument must be a
            % vector with a number of components that is an integral multiple of 4,
            % i.e. in our case it must have 4 components:
            Screen('DrawTexture', win, gratingtex, [], [], stimulus, [], [], [], [], rotateMode, [phase, freq, actualAmplitude, 0]);
            
            % Show it at next retrace:
            if k==1
                vbl = Screen('Flip', win, vbl + randomNumber);
            else
                vbl = Screen('Flip', win, vbl + flashDuration);
            end
            if useArduino, writeDigitalPin(a, 'D46', 1);  end % send arduino signal 0 (empty screen)
            data(currentLine,1)=datenum(datetime('now','Format','MM/dd/yy HH:mm:ss.SSSSSS'));
            data(currentLine,3)= flashStatus; % flipped compared to the previous one
            if write, dlmwrite(filename,data,'precision','%.8f'); end
            flashStatus = 1-flashStatus;
            figure(1); clf; plot(data(data(:,2)>-2,2),'.'); ylim(ylim+[-1 1]*30); PlotHVLines(sum(data(1:currentLine-1,2)>-2),'v','r');
            xlabel('trial ID'); ylabel('orientation (degrees)'); drawnow
            currentLine = currentLine+1;
        end
    end
end
% We're done. Close the window. This will also release all other ressources:
sca;

beep

% Bye bye!
return;
