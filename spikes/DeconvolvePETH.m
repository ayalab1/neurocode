function [deconvolved,t] = DeconvolvePETH(signal,events,varargin)

% [DeconvolvePETH] - [Compute a deconvolved version of PETH which removes the
% smoothing effects of the events' autocorrelogram]
%
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% default values:
durations = [-1 1];
nBins = [100];

argsToPassOn = {};
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PETH">PETH</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'durations'
            durations = varargin{i+1};
            dt = durations(2);
            if ~isvector(dt) || length(dt) ~= 1
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'duration'
            durations = varargin{i+1};
            dt = durations(2);
            if ~isvector(dt) || length(dt) ~= 1
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'nbins'
            nBins = varargin{i+1};
            if ~isvector(nBins) || length(nBins) ~= 1
                error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'smooth'
            smooth = varargin{i+1};
            if ~isvector(smooth) || length(smooth) ~= 1
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        otherwise
            argsToPassOn{end+1} = varargin{i}; argsToPassOn{end+1} = varargin{i+1};
    end
end

[~,t] = PETH(signal,events,'durations',durations,'nBins',nBins,argsToPassOn{:});
[autoResponse,~] =  PETH(events,events,'durations',durations*2,'nBins',nBins*2-1,argsToPassOn{:});
[rawResponse,~] = PETH(signal,events,'durations',durations*2,'nBins',nBins*2-1,argsToPassOn{:}); % make sure the binsize stays the same while you double the duration

autocorrelogram = sum(autoResponse); rawPeth = sum(rawResponse); 
const = mean(rawPeth); rawPeth = rawPeth - const; % remove the mean, because not all firing needs to be explained by the stimulus
T0 = toeplitz([autocorrelogram(:); zeros(numel(rawPeth)-numel(autocorrelogram), 1)], [autocorrelogram(1), zeros(1, length(autocorrelogram)-1)]);
T = T0(nBins:end,1:nBins);
deconvolved = T \ rawPeth(nBins/2+0:nBins/2*3-1)' + const/length(events); % add the baseline to the final PETH


%% Old code:
% function deconvolved = DeconvolvePETH(rawPeth,autocorrelogram)
% 
% Note: "rawPETH" should be of double the duration of the autocorrelogram for this to work/

% Example:
% [autoResponse,t]=  PETH(events,events,'durations',durations,'nBins',nBins);
% [rawResponse,~]=  PETH(signal,events,,'durations',durations*2,'nBins',nBins*2-1); % make sure the binsize stays the same while you double the duration
% deconvolved = DeconvolvePETH(sum(rawResponse),sum(autoResponse));
% plot(t,deconvolved);
% 
% T = toeplitz([autocorrelogram(:); zeros(numel(rawPeth)-numel(autocorrelogram), 1)], [autocorrelogram(1), zeros(1, length(autocorrelogram)-1)]);
% deconvolved = T \ rawPeth(:);