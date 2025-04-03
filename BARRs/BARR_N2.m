function BARR_N2(varargin)
% barNeuroScope2: make single event file for barrages (HSEn2) in specified basepath
%
% Syntax:
%   barNeuroScope2('basepath', basepath)
%
% Inputs:
%   basepath - (optional) string, path to the directory containing the input
%              files and where the output file will be saved. Default: current
%              working directory (pwd).
%   clean - (optional) logical, whether or not to use clean trial
%
% Outputs:
%   Saves a structure called HSEn2 to a file called [basename '.HSEn2.events.mat']
%   in the specified basepath. The HSEn2 structure has two fields:
%   - timestamps: an array of timestamps for barrage events.
%   - peaks: an array of peak amplitudes for barrage events.
%
% Notes:
%   - The input files 'HSEfutEVT.mat' and 'HSE.mat' must be located in the
%     subdirectory 'Barrage_Files' within the specified basepath.
%   - The function uses the variables evtSave, HSE, HSE.keep, and HSE.NREM,
%     which are expected to be defined in the input files.
%   - The timestamps and peaks arrays are created from the evtSave variable,
%     using only the rows corresponding to the HSE.keep(HSE.NREM) indices.

p=inputParser;
addParameter(p,'basepath',pwd); 
addParameter(p,'clean',0);
parse(p,varargin{:});
basepath = p.Results.basepath;
clean = p.Results.clean;

%% Create Neuroscope2 file for barrage events
basename = basenameFromBasepath(basepath);

savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
if clean
    load([savePath 'clean.HSE.mat']);
else
    load([savePath 'HSE.mat']);
end

HSEn2.timestamps = HSE.timestamps(HSE.keep(HSE.NREM),:);
HSEn2.peaks = HSE.peaks(HSE.keep(HSE.NREM),:);
HSEn2.duration = HSEn2.timestamps(:,2) - HSEn2.timestamps(:,1);

save(fullfile(basepath,[basename '.HSEn2.events.mat']), 'HSEn2');
end