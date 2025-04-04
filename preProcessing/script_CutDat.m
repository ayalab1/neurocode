% script_CutDat
% Use this script if you want to remove the end of your dat files (ephys and other input files).
% Just set the number of channels and the cutoff time (max time) (in seconds)
% in the lines below and run the script!
% The original files will remain but renamed to append _original. Feel free
% to remove them afterwards.
%
% Copyright (C) 2021-2023 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


nChannels = ?; % set this manually
maxTime =; % set manually: this should be in seconds

% rename files: do this only once as doing it more than once would remove the original files!!!
if exist('amplifier_original.dat', 'file')
    error('amplifier_original.dat file already exists. Don''t run script again, it might overwrite your original files!!!')
else
    movefile('auxiliary.dat', 'auxiliary_original.dat');
    movefile('digitalin.dat', 'digitalin_original.dat');
    movefile('analogin.dat', 'analogin_original.dat');
    movefile('time.dat', 'time_original.dat');
    movefile('amplifier.dat', 'amplifier_original.dat');
end

file = memmapfile('amplifier_original.dat', 'Format', 'int16', 'Writable', false);
data = reshape(file.data, nChannels, []);
originalSize = size(data, 2);
data = data(:, 1:ceil(maxTime*20000));
f = fopen('amplifier.dat', 'w'); % new file should be test.dat
fwrite(f, data(:), 'int16');
fclose(f);

file = memmapfile('auxiliary_original.dat', 'Format', 'uint16', 'Writable', false);
nChannels = length(file.data) / originalSize;
data = reshape(file.data, nChannels, []);
data = data(:, 1:ceil(maxTime*20000));
f = fopen('auxiliary.dat', 'w'); % new file should be test.dat
fwrite(f, data(:), 'uint16');
fclose(f);

file = memmapfile('digitalin_original.dat', 'Format', 'uint16', 'Writable', false);
nChannels = length(file.data) / originalSize;
data = reshape(file.data, nChannels, []);
data = data(:, 1:ceil(maxTime*20000));
f = fopen('digitalin.dat', 'w'); % new file should be test.dat
fwrite(f, data(:), 'uint16');
fclose(f);

file = memmapfile('time_original.dat', 'Format', 'int32', 'Writable', false);
nChannels = length(file.data) / originalSize;
data = reshape(file.data, nChannels, []);
data = data(:, 1:ceil(maxTime*20000));
f = fopen('time.dat', 'w'); % new file should be test.dat
fwrite(f, data(:), 'int32');
fclose(f);

file = memmapfile('analogin_original.dat', 'Format', 'uint16', 'Writable', false);
nChannels = length(file.data) / originalSize;
data = reshape(file.data, nChannels, []);
data = data(:, 1:ceil(maxTime*20000));
f = fopen('analogin.dat', 'w'); % new file should be test.dat
fwrite(f, data(:), 'uint16');
fclose(f);
