InitializePsychSound(1); nrchannels = 1;
pahandle = PsychPortAudio('Open', [], [], 0, round(1/mode(diff(sound(:,1)))), 1);
PsychPortAudio('FillBuffer', pahandle, soundDownward(:,2)');
t1 = PsychPortAudio('Start', pahandle, 1, 0, 1);

%% 

InitializePsychSound(1); nrchannels = 1;
pahandle = PsychPortAudio('Open', [], [], 0, round(1/mode(diff(sound(:,1)))), 1);
PsychPortAudio('FillBuffer', pahandle, soundUpward(:,2)');
t1 = PsychPortAudio('Start', pahandle, 1, 0, 1);
