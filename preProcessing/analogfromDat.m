function analogfromDat(basepath,varargin)
%
%
%   [perform lowpass (2 X output Fs) sinc filter on wideband data
%   (specifically converting analog into downsampled signal similar to lfp
%   from DAT extraction). Specific use case is extracting sinusoids from
%   opto stim for CFCPhaseAmp analysis. subsample the filtered data and save as a new flat binary
%   Must have analogin.dat file in the basepath that is concatenated
%   analogin.dat
%
%   [note that sincFilter was altered to accomodate GPU filtering]
%
% INPUTS
% [basePath]    [path where the recording files are located
%               Assumes presence of the following files:
%                   basePath/analogin.dat  
%
%                   (optional parameters files)
%                   basePath/baseName.xml
%                   basePath/baseName.sessionInfo.mat
%
%               If basePath not specified, tries the current working
%               directory]
%   (options)
%    =========================================================================
%     Properties    Values
%    ------------------------------------------------------------------------
%       ['outFs']       [(default: 1250) downsampled frequency of the .lfp
%                       output file. if no user input and not specified in
%                       the xml, use default]
%       ['lopass']      [(default: 450) low pass filter frequency]
%       ['useGPU']      [(default: false) whether or not to use GPU to speed
%                       processing (might not want to if limited GPU)]
%
%
%  OUTPUT
%   [Creates file:   basePath/baseName.analogin]
%
%   [If no sessionInfo.mat file previously exists, creates one with
%   the information from the .xml file, with the .analogin sampling frequency
%   and the lowpass filter used]
%
% SEE ALSO
%
% Dependency: iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
%
%[HLarsson] [2023-2024]
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.




%% Import

import iosr.dsp.*

%% Input handling

if ~exist('basepath','var')
    basepath = pwd;
end

p = inputParser;

addParameter(p,'outFs',1250,@isnumeric);
addParameter(p,'lopass',450,@isnumeric);
addParameter(p,'useGPU',false,@islogical);



parse(p,varargin{:})
outFs = p.Results.outFs;
lopass = p.Results.lopass;
useGPU = p.Results.useGPU;


session = getSession('basepath',basepath);
basename = session.general.name;


datFile = 'analogin.dat';


% if no gpuDevice found dont try to use
if useGPU && gpuDeviceCount<1
    warning('No GPU device found, continuing without..')
    useGPU = false;
end

if useGPU
    g = gpuDevice(1);
end

sizeInBytes = 2; %

%% housekeeping

%Check the dat
fInfo = checkFile('basepath',basepath,'filename',datFile,'searchSubdirs',false');
fdat = fInfo.name;

%Get the metadata - update approach to avoid descripenancies if session
%file has different sampling rate (would theoritically be fine but still)

[amplifier_channels, notes, aux_input_channels, spike_triggers,...         
        board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
        read_Intan_RHD2000_file_bz;
inFs = frequency_parameters.amplifier_sample_rate;

[a, aa] = size(aux_input_channels);
nbChan = aa;

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
flfp = fullfile(basepath,[basename,'.analogin']);


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
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunksize))-1;

%% GET LFP FROM DAT

if exist([basepath '\' basename '.analogin'],'file') 
    fprintf('.analogin file already exists \n')
    return
end

fidI = fopen(fdat, 'r');
fprintf('Extraction of analog signal has begun \n')
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
        try
            dat = reshape(dat,[nbChan (chunksize+2*ntbuff)]);
        catch
            % One possible issue is that the network bugged and so the file was dropped. This can be fixed by reloading the file:
            % === SOLUTION ===
            fidI = fopen(fdat, 'r');
            fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');
            dat = fread(fidI,nbChan*(chunksize+2*ntbuff),'int16');
            % === END OF SOLUTION === % if this executes fine, hit "dbcont"
            try
                dat = reshape(dat,[nbChan (chunksize+2*ntbuff)]);
            catch
                warning('This should be fixed.. tell Raly! [or repeat the solution just above and see if that works. If no errors, hit dbcont!]');
                keyboard;
            end
        end
    else
        dat = fread(fidI,nbChan*(chunksize+ntbuff),'int16');
        try
            dat = reshape(dat,[nbChan (chunksize+ntbuff)]);
        catch
            warning('This should be fixed.. tell Raly!');
            keyboard;
        end
    end
    
    
    DATA = nan(size(dat,1),chunksize/sampleRatio);
    for ii = 1:size(dat,1)
        
        d = double(dat(ii,:));
        if useGPU
            d = gpuArray(d);
            tmp = gpuArray(zeros(size(d)));
        end
        
        tmp = iosr.dsp.sincFilter(d,ratio);
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
    try
        dat = reshape(dat,[nbChan (remainder+ntbuff)]);
    catch
        warning('Check the number of channels in the xml match what you recorded. If not the problem, tell Raly!');
        keyboard;
    end
    
    DATA = nan(size(dat,1),floor(remainder/sampleRatio));
    for ii = 1:size(dat,1)
        d = double(dat(ii,:));
        if useGPU
            d = gpuArray(d);
        end
        
        tmp = iosr.dsp.sincFilter(d,ratio);
        
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

if useGPU
    reset(g);
    gpuDevice([]);
end


end