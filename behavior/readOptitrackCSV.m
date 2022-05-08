function [tracking] = readOptitrackCSV(filename,varargin)

% INPUTS
%    fbasename   -basename of the recording (this is only used in generating
%                 the .mat output file)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'syncDatFile' name of binary file where sync signal is stored
%                   (default = filebasename_digitalin.dat)
%     'syncSampFq'  sampling freqeuncy of the sync signal (default= 20kHz)
%     'syncChan'    sync channel to read (default = 1)
%     'syncNbCh'    number of channels in the sync file (default = 1)
%     'posSampFq'   sampling frequency of the final pos file (after
%                   interpolation) in Hz
%     'columnOrder' order of output variales in .csv file: frame, time, rx,
%                   ry, rz, rw, x, y, z, error. Vector with postion of each
%                   variable in the order specified here
%     'basepath'    absolute path to folder with .csv optitrack files (default = pwd) 
%
% OUTPUTS
%   tracking - strcutre with data read from csv file and timestamps

%% Parse Inputs
p = inputParser;

addParameter(p,'syncDatFile','digitalin.dat',@ischar)
addParameter(p,'syncSampFq',20000,@isnumeric)
addParameter(p,'syncChan',1,@isnumeric)
addParameter(p,'syncNbCh',1,@isnumeric)
addParameter(p,'posSampFq',120,@isnumeric)
addParameter(p,'columnOrder',1:9,@isnumeric)
addParameter(p,'basepath',pwd,@ischar)

parse(p,varargin{:});
syncDatFile = p.Results.syncDatFile;
syncSampFq = p.Results.syncSampFq;
syncChan = p.Results.syncChan;
syncNbCh = p.Results.syncNbCh;
posSampFq = p.Results.posSampFq;
columnOrder = p.Results.columnOrder;
basepath = p.Results.basepath;

% Check tracking file
checkFile('basepath',basepath,'filename',filename,'fileType','.csv');

%% Import and correct data

dat = importdata([basepath, filesep, filename]);
% dat = bz_scrubTracking(dat); % modified by AFR
pos = dat.data; % all tracking variables   

% read optitrack sync channel
fid = fopen([basepath, filesep, syncDatFile]);
dig = fread(fid,[syncNbCh inf],'int16=>int16');  % default type for Intan digitalin
dig = dig(syncChan,:);
t = (0:length(dig)-1)'/syncSampFq; % time vector in sec


% get frame timing in digital input timestamps
dPos = find(diff(dig)==1);
dNeg = find(diff(dig)==-1);

if length(dPos) == length(dNeg)+1
    dPos = dPos(1:end-1);
elseif length(dNeg) == length(dPos)+1
    dNeg = dNeg(1:end-1);
elseif abs(length(dNeg)-length(dPos)) > 1
    warning('some problem with frames');
    keyboard
end
% Frame timing is the middle of shuter opening
frameT = (t(dPos)+t(dNeg))/2;


% The system sometimes (rarely) keeps on recording a few frames after software stopped
% recording. So we skip the last frames of the TTL
if length(frameT)<size(pos,1)
    warning('Too many video frames!'); % maybe because intan was stopped before Optitrack
    %keyboard
    pos(pos==-1) = NaN;
    pos = pos(1:length(frameT),:); % ???
elseif length(frameT)>size(pos,1)
    frameT = frameT(1:size(pos,1));
end

% We now interpolate the data at 120 Hz (or sampling fqcy specified in arguments)
recDuration = length(dig)/syncSampFq;

timestamps = (0:1/posSampFq:recDuration-1/posSampFq)';
pos(pos==-1) = NaN;
newPos = interp1(frameT,pos,timestamps);
newPos = newPos(:,columnOrder);

tracking.timestamps = timestamps;
tracking.data = newPos;

end
