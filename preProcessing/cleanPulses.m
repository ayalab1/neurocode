function [data] = cleanPulses(ts, varargin)
%   
%  [[data] = cleanPulses(ts, varargin)]
%
% [Remove artifacts form square-shappped optogetic pulses form dat file]
%
% INPUTS
%
%  [ts]        [List of artifacts that will be removed from data (in
%               seconds)]
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     ['filename']  [filename to clean]
%     ['basepath']  [basepath to clean]
%     ['correctDC'] [Logical variable to indicate if DC is corrected, default
%                   false]
%     ['ch']        [List of channels to clean pulses, default all]
%     ['winArt']    [window for artefact removal, in seconds, default
%                   0.0005s]
%     ['winDC']     [window for DC removal, in seconds, default 0.005s]
%
%  OUTPUTS
%  [data]   [cleaned dat file]
%
%  SEE ALSO
%
% [Manu Valero, Antonio FR] [2021-2022]
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default parameters
filename = split(pwd,filesep); filename = filename{end};

% Parse options
p = inputParser;
addParameter(p,'filename',filename,@isstr);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'correctDC',false, @islogical);
addParameter(p,'ch','all');
addParameter(p,'winArt',0.0005,@isnumeric);
addParameter(p,'winDC',0.005,@isnumeric);

parse(p, varargin{:});

filename = p.Results.filename;
basepath = p.Results.basepath;
correctDC = p.Results.correctDC;
ch = p.Results.ch;
winArt = p.Results.winArt;
winDC = p.Results.winDC;

try [session] = getSession('basepath',basepath);
    fs = session.extracellular.sr;
    nChannels = session.extracellular.nChannels;
catch
    warning('SessionInfo file not found.');
end

if ischar('ch') && strcmpi(ch, 'all')
    ch = 1:nChannels;
else
    ch = ch + 1;
end

ts = int32(ts * fs);
winArt = winArt * fs;
winDC = winDC * fs;

% disp('Creating .dat back up...');
% copyfile(strcat(filename,'.dat'),'copy_bin.dat'); % create .dat back up
m = memmapfile(fullfile(basepath,strcat(filename,'.dat')),'Format','int16','Writable', true);
data=reshape(m.Data,nChannels,[]);

for hh = ch
    fprintf('Channel #%i of %i...\n',hh,length(ch))
    for ii = 1:size(ts,2) % tsses
        % remove dc
        if correctDC
            data(hh,ts(1,ii):ts(2,ii)) = data(hh,ts(1,ii):ts(2,ii)) -...
                median(data(hh,ts(1,ii):ts(1,ii)+winDC)) + ...
                median(data(hh,ts(1,ii)-winDC:ts(1,ii)));
        end
        % remove artifacts
        for jj = 1:size(ts,1)
            art = int32([ts(jj,ii) - winArt ts(jj,ii) + winArt]);
            data(hh,art(1):art(2))=int16(interp1(double(art),...
                double(data(hh,art)),double(art(1):art(2))));
        end
    end
end
disp('Cleaning data on disk...');
m.Data=data(:);

end