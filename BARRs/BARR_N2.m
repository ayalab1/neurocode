function BARR_N2(basepath)
% barNeuroScope2: make single event file for barrages (HSEn2) in specified basepath
%
% Inputs:
%   basepath - (optional) string, path to the directory containing the input
%              files and where the output file will be saved. Default: current
%              working directory (pwd).
%
% Outputs:
%   Saves a structure called HSEn2 to a file called [basename '.HSEn2.events.mat']
%   in the specified basepath. The HSEn2 structure has three fields, all in seconds:
%   - timestamps: an array of timestamps for barrage events.
%   - peaks: an array of peak amplitudes for barrage events.
%   - duration: array of duration for barrage events.
%
% Notes:
%   - 'HSE.mat' must be located in the
%     subdirectory 'Barrage_Files' within the specified basepath.

if nargin<1
    basepath = pwd;
end

%% Create Neuroscope2 file for barrage events
basename = basenameFromBasepath(basepath);

savePath = strcat(basepath, filesep, 'Barrage_Files', filesep, basename, '.');
load([savePath 'HSE.mat']);

HSEn2.timestamps = HSE.timestamps(HSE.keep(HSE.NREM),:);
HSEn2.peaks = HSE.peaks(HSE.keep(HSE.NREM),:);
HSEn2.duration = HSEn2.timestamps(:,2) - HSEn2.timestamps(:,1);

save(fullfile(basepath,[basename '.HSEn2.events.mat']), 'HSEn2');
end