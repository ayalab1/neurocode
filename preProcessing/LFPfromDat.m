function LFPfromDat(basepath,varargin)
% perform lowpass (2 X output Fs) sinc filter on wideband data
% subsample the filtered data and save as a new flat binary
% basename must have basename.dat and basename.xml
% basepath is the full path for basename.dat
%
% note that sincFilter was altered to accomodate GPU filtering
%
%INPUTS
%   basePath    path where the recording files are located
%               where basePath is a folder of the form: 
%                   whateverPath/baseName/
%
%               Assumes presence of the following files:
%                   basePath/baseName.dat
%                   -or-
%                   basePath/amplifier.dat
%
%                   (optional parameters files)
%                   basePath/baseName.xml
%                   basePath/baseName.sessionInfo.mat
%
%               If basePath not specified, tries the current working directory.
%
%   (options)
%       'datFile'       specify while file you'd like to compute lfp from
%       'outFs'         (default: 1250) downsampled frequency of the .lfp
%                       output file. if no user input and not specified in
%                       the xml, use default
%       'lopass'        (default: 450) low pass filter frequency 
%       'useGPU'        (default: false) whether or not to use GPU to speed
%                       processing (might not want to if limited GPU)
%
%
%OUTPUT
%   Creates file:   basePath/baseName.lfp
%
%   If no sessionInfo.mat file previously exists, creates one with 
%   the information from the .xml file, with the .lfp sampling frequency 
%   and the lowpass filter used.
%
%
%Dependency: iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
%
%SMckenzie, BWatson, DLevenstein 2018
%
%kathrynmcclain 2020

%% Import

import iosr.dsp.*

%% Input handling

if ~exist('basepath','var')
    basepath = pwd;
end

p = inputParser;
addParameter(p,'datFile',[],@isstr);
addParameter(p,'outFs',[],@isnumeric);
addParameter(p,'lopass',450,@isnumeric);
addParameter(p,'useGPU',false,@islogical);
addParameter(p,'inFs',[],@isnumeric);

parse(p,varargin{:})
datFile = p.Results.datFile;
outFs = p.Results.outFs;
lopass = p.Results.lopass;
useGPU = p.Results.useGPU;
inFs = p.Results.inFs;

session = getSession('basepath',basepath);
basename = session.general.name;

if isempty(datFile)
    datFile = [basename,'.dat'];
elseif ~strcmp(datFile(end-3:end),'.dat')
    datFile = [datFile,'.dat'];
end

% if no gpuDevice found dont try to use
if useGPU && gpuDeviceCount<1
    warning('No GPU device found, continuing without..')
    useGPU = false;
end

sizeInBytes = 2; %

%% housekeeping

%Check the dat
fInfo = checkFile('basepath',basepath,'filename',datFile);
fdat = [fInfo.folder, filesep, fInfo.name];

%Get the metadata
if isempty(inFs)
    inFs = session.extracellular.sr;
end

nbChan = session.extracellular.nChannels;

%set output sampling rate from xml, user input
if isempty(outFs)          %If user input - priority (keep from above)
    outFs = session.extracellular.srLfp;
end

if lopass> outFs/2
    warning('low pass cutoff beyond Nyquist')
end
 
ratio = lopass/(inFs/2) ;
sampleRatio = (inFs/outFs);

%output file
flfp = fullfile(basepath,[basename,'.lfp']);

%% Set Chunk and buffer size at even multiple of sampleRatio
chunksize = 1e5; % depends on the system... could be bigger I guess
if mod(chunksize,sampleRatio)~=0
    chunksize = chunksize + sampleRatio-mod(chunksize,sampleRatio);
end

%ntbuff should be even multiple of sampleRatio
ntbuff = 525;  %default filter size in iosr toolbox
if mod(ntbuff,sampleRatio)~=0
    ntbuff = ntbuff + sampleRatio-mod(ntbuff,sampleRatio);
end

nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunksize));

%% GET LFP FROM DAT

if exist([basepath '\' basename '.lfp'])
    fprintf('LFP file already exists \n')
else
fidI = fopen(fdat, 'r');
fprintf('Extraction of LFP begun \n')
fidout = fopen(flfp, 'a');

for ibatch = 1:nbChunks

    if mod(ibatch,10)==0
        if ibatch~=10
            fprintf(repmat('\b',[1 length([num2str(round(100*(ibatch-10)/nbChunks)), ' percent complete'])]))
        end
        fprintf('%d percent complete', round(100*ibatch/nbChunks));
    end
    
    if ibatch>1
        fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');
        dat = fread(fidI,nbChan*(chunksize+2*ntbuff),'int16');
        dat = reshape(dat,[nbChan (chunksize+2*ntbuff)]);
    else
        dat = fread(fidI,nbChan*(chunksize+ntbuff),'int16');
        dat = reshape(dat,[nbChan (chunksize+ntbuff)]);
    end
    
    
    DATA = nan(size(dat,1),chunksize/sampleRatio);
    for ii = 1:size(dat,1)
        
        d = double(dat(ii,:));
        if useGPU
            d = gpuArray(d);
            tmp = gpuArray(zeros(size(d)));
        end
        
        tmp=  iosr.dsp.sincFilter(d,ratio);
        if useGPU
            if ibatch==1
                DATA(ii,:) = gather_try(int16(real( tmp(sampleRatio:sampleRatio:end-ntbuff))));
            else
                DATA(ii,:) = gather_try(int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff))));
            end
            
        else
            if ibatch==1
                DATA(ii,:) = int16(real( tmp(sampleRatio:sampleRatio:end-ntbuff)));
            else
                DATA(ii,:) = int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff)));
            end
            
        end
        
    end
    
    fwrite(fidout,DATA(:),'int16'); 
end



remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunksize;
if ~isempty(remainder)
    fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');
    dat = fread(fidI,nbChan*(remainder+ntbuff),'int16');
    dat = reshape(dat,[nbChan (remainder+ntbuff)]);
 
    DATA = nan(size(dat,1),floor(remainder/sampleRatio));
    for ii = 1:size(dat,1)
        d = double(dat(ii,:));
        if useGPU
            d = gpuArray(d);
        end
        
        tmp=  iosr.dsp.sincFilter(d,ratio);
        
        if useGPU
            
            DATA(ii,:) = gather_try(int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end))));
        else
            DATA(ii,:) = int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end)));
        end
    end
    
    fwrite(fidout,DATA(:),'int16');
end

fclose(fidI);
fclose(fidout);

end
disp('lfp file created')
end

