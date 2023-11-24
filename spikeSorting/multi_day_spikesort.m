function multi_day_spikesort(animal, sessions, shanks)
% multi_day_spikesort: creates individual dat file over sessions for shank/s
%
%
% Example
%
% animal = 'Z:\Data\HMC1';
% sessions = {'day8', 'day9', 'day10'};
% shanks = [1,2,3,4];
% multi_day_spikesort(animal,sessions,shanks)

batch_size = 20000 * 10; % samples

% load session metadata
load(fullfile(animal, sessions{1}, [sessions{1}, '.session.mat']), 'session')

% get channels we want to use/concatenate
chans = [];
for shank = shanks'
    chans = [chans, session.extracellular.electrodeGroups.channels{shank}];
end
chans = sort(chans);

% memmap dat files
for session_i = 1:length(sessions)
    current_dat_path = fullfile(animal, sessions{session_i}, [sessions{session_i}, '.dat']);
    d = dir(current_dat_path);
    nSamp = d.bytes / 2 / session.extracellular.nChannels;
    mmf{session_i} = memmapfile(current_dat_path, 'Format', {'int16', [session.extracellular.nChannels, nSamp], 'x'});
end

% write new dat file within new folder
mkdir(fullfile(animal, [sessions{:}]))

f = fopen(fullfile(animal, [sessions{:}], [[sessions{:}], '.dat']), 'a');

% anonymous function to flatten signal for saving
% flatten = @(x) x(:);

% iterate over sessions
for session_i = 1:length(sessions)
    disp(['day ', num2str(session_i), ' of ', num2str(length(sessions))])

    n_samples = size(mmf{session_i}.Data.x, 2);

    batches = round(linspace(1, size(mmf{session_i}.Data.x, 2), round(n_samples/batch_size)));

    WaitMessage = parfor_wait(length(batches)-1);

    % iterate over batches
    for batch_i = 1:length(batches) - 1
        WaitMessage.Send;
        % fwrite(f, flatten(mmf{session_i}.Data.x(chans, batches(batch_i):batches(batch_i+1))), 'int16');
        fwrite(f, mmf{session_i}.Data.x(chans, batches(batch_i):batches(batch_i+1)), 'int16');

    end
    WaitMessage.Destroy;
end
fclose(f);

end