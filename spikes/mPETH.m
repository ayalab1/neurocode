function [m,t] = mPETH(data, event, varargin)

% [m,t] = mPETH(data,events,<options>).
% returns the mean PETH (see <a href="matlab:help PETH">PETH</a>)
% Copyright (C) 2018 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% default values
duration = [-1 1];
pictureoptions = {};
smooth = 0;
nBins = 100;

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PETH">PETH</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'durations',
            duration = varargin{i+1};
            dt = duration(2);
            if ~isvector(dt) || length(dt) ~= 1,
                error('Incorrect value for property ''dt'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'duration',
            duration = varargin{i+1};
            dt = duration(2);
            if ~isvector(dt) || length(dt) ~= 1,
                error('Incorrect value for property ''dt'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'nbins',
            nBins = varargin{i+1};
            if ~isvector(nBins) || length(nBins) ~= 1,
                error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'show',
            show = varargin{i+1};
            if ~isastring(show,'on','off'),
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'mode',
            mode = varargin{i+1};
            if ~isastring(mode,'l','c'),
                error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'name',
            namestring = varargin{i+1};
            if ~isastring(namestring),
                error('Incorrect value for property ''name'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'smooth',
            smooth = varargin{i+1};
            if ~isvector(smooth) || length(smooth) ~= 1,
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'linetype',
            linetype = varargin{i+1};
            if ~isastring(linetype),
                error('Incorrect value for property ''linetype'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        otherwise,
            pictureoptions{end+1} = varargin{i};%pictureoptions{end+1} = varargin{i+1};
    end
end

if isempty(event)
    m = nan(1,nBins);
    t = linspace(duration(1),duration(2),nBins);
    return
end
if size(data,2)==1 && exist('fastPETH')
    peths = fastPETH(data,event,duration);
    j = size(event,1)*ones(size(peths));
else
    [peths, j] = Sync(data, event, 'durations', duration);
end
[m,~,t] = SyncHist(peths,j, 'nBins', nBins, 'smooth', smooth, 'mode', 'mean', 'durations', duration);
m = m(:)';
if isempty(m)
    m = zeros(1,nBins);
end



