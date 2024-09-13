function cut_dat(nChannels, maxTime, varargin)
% cut_dat: Trims the end of specified .dat files (ephys and other input files).
%
% This function removes excess data from the end of various .dat files 
% commonly used in electrophysiology recordings (e.g., amplifier, auxiliary, 
% digitalin). It retains the original files by renaming them with an "_original" 
% suffix, allowing you to manually delete them afterward if desired.
%
% Inputs:
%   - nChannels: The number of channels in the data files (e.g., 128).
%   - maxTime: The maximum time (in seconds) to keep in the trimmed file.
%   - varargin: Optional parameters, including:
%       - 'sample_rate': The sample rate of the recordings (default is 20,000 Hz).
%
% Usage:
%   Change to the directory containing the .dat files and run:
%     cut_dat(nChannels, maxTime)
%   Example:
%     cd('path_to_data')
%     cut_dat(128, 14586.933)
%
% Functionality:
%   - Renames the original files by appending "_original" to the filename. 
%     This only happens once to prevent overwriting.
%   - Processes the .dat files in batches to avoid reading the entire file into memory.
%   - Maintains the correct data type for each file (e.g., int16, uint16).
%
% Notes:
%   - This function is based on "script_CutDat.m" but improves performance 
%     by using batching, proper data types, and handling various file types.
%   - Ensure the renamed "_original" files are not accidentally deleted until 
%     you're certain the trimming has been successful.
%
% Errors:
%   - If the script is run multiple times, it will throw an error if 
%     '_original.dat' files already exist, preventing accidental overwriting.
%
% Ryan
% 

% parse inputs
p = inputParser;
addParameter(p, 'sample_rate', 20000, isnumeric(x));

parse(p, varargin{:})
sample_rate = p.Results.sample_rate;


files_table = table();
files_table.files = {'amplifier', 'auxiliary', 'digitalin', 'digitalout', 'analogin', 'time', 'supply'}';
files_table.data_type = {'int16', 'uint16', 'uint16', 'uint16', 'uint16', 'int32', 'uint16'}';
disp(files_table)

% rename files: do this only once as doing it more than once would remove the original files!
if exist('amplifier_original.dat', 'file')
    error('amplifier_original.dat file already exists. Don''t run script again, it might overwrite your original files!!!')
else
    % iter over each file type and make a copy
    for i = 1:length(files_table.files)
        if exist([files_table.files{i}, '.dat'], 'file')
            movefile([files_table.files{i}, '.dat'], [files_table.files{i}, '_original.dat'])
        end
    end
end

for i = 1:length(files_table.files)

    if ~exist([files_table.files{i}, '_original.dat'], 'file')
        continue
    end
    disp(files_table.files{i})

    % first get n samples from amplifier and use the originalSize for the
    % remaining files
    if contains(files_table.files{i}, 'amplifier')
        d = dir([files_table.files{i}, '_original.dat']);
        originalSize = d.bytes / 2 / nChannels;

    else
        d = dir([files_table.files{i}, '_original.dat']);
        nChannels = d.bytes / 2 / originalSize;
    end

    disp('writing file...')
    DAT_extractChannel([files_table.files{i}, '_original.dat'], [files_table.files{i}, '.dat'], ...
        nChannels, ...
        1:nChannels, ...
        0, ...
        ceil(maxTime*sample_rate), ...
        files_table.data_type{i}, ...
        0, ...
        0)
end
end

function DAT_extractChannel(fname, fname_new, numchannel, chselect, read_start, read_until, precision, b_skip, append_mode)
% eeg is output, fb is size of file in bytes
% Reads multi-channel recording file to a matrix
% last argument is optional (if omitted, it will read all the
% channels.
%
% From the Buzsaki lab (obtained 4/5/2010).
% revised by zifang zhao @ 2014-5-1 increased 2 input to
% control the range of file reading

if nargin < 7 %precision and skip
    precision = 'int16';
end
if nargin < 8 %skip
    b_skip = 0;
end
if nargin < 9
    append_mode = 0;
end

if isempty(precision)
    precision = 'int16';
end
if isempty(b_skip)
    b_skip = 0;
end
if append_mode == 0
    mode = 'w+';
else
    mode = 'a+';
end
fileinfo = dir(fname);
if nargin == 3
    datafile = fopen(fname, 'r');
    eeg = fread(datafile, [numchannel, inf], precision);
    fclose(datafile);
    eeg = eeg';
    return
end

fh_write = fopen(fname_new, mode);

fb = fileinfo(1).bytes;
numel_all = floor(fb/2/numchannel);
if nargin >= 3
    % the real buffer will be buffersize * numch * 2 bytes
    % (int16 = 2bytes)
    if nargin < 4
        read_until = numel_all;
    end
    buffersize = 4096;
    % get file size, and calculate the number of samples per channel
    if read_start < 0
        read_start = read_start + numel_all - 1;
        if read_until == 0
            read_until = numel_all;
        end
    end

    read_start(read_start < 0) = 0;
    read_until = read_until + 1;
    read_until(read_until > numel_all) = numel_all;
    read_start_byte = read_start * 2 * numchannel;
    read_until_byte = read_until * 2 * numchannel;

    %% original method
    numel1 = 0;
    datafile = fopen(fname, 'r');
    state = fseek(datafile, read_start_byte, 'bof');
    while ftell(datafile) < read_until_byte && ~feof(datafile) && state == 0
        len_left = read_until_byte - ftell(datafile);
        if len_left >= buffersize * numchannel * 2
            [data, count] = fread(datafile, [numchannel, buffersize], precision, b_skip); %can be improved,vectorize,arrayfun,multi-threading, zifangzhao@4.24
        else
            [data, count] = fread(datafile, [numchannel, ceil(len_left/numchannel/2)], precision, b_skip); %can be improved,vectorize,arrayfun,multi-threading, zifangzhao@4.24
        end
        numelm = (count / numchannel);
        if numelm > 0

            fwrite(fh_write, data(chselect, :), precision);
            numel1 = numel1 + numelm;
        end
    end
end
fclose(datafile);
fclose(fh_write);
end