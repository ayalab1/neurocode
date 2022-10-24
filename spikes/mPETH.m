function [m,t] = mPETH(data, event, varargin)
%
% [[m,t] = mPETH(data,events,<options>)]
% [returns the mean PETH (see <a href="matlab:help PETH">PETH</a>)]
%
%
%
%  INPUTS
%
%    [data]        [data to find events for]
%    [events]      [events for PETH]
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-1 1])
%     'nBins'       number of time bins around the events (default = 101)
%     'mode'        whether the sample data is linear ('l') or circular ('c')
%                   (for example, in the case 'samples' is the phase of an
%                   oscillation).
%     'show'        display the mean PETH (default = 'on' when no outputs are
%                   requested and 'off' otherwise)
%     'smooth'      standard deviation for Gaussian kernel (default = 1 bin)
%                   applied to the mean peri-event activity 'm' (note, no
%                   smoothing will be applied to the output 'matrix')
%     'title'       if the results are displayed ('show' = 'on'), specify a
%                   desired title (default is deduced by variable names)
%     <plot options> any other property (and all the properties that follow)
%                   will be passed down to "plot" (e.g. 'r', 'linewidth', etc)
%                   Because all the following inputs are passed down to "plot",
%                   make sure you put these properties last.
%    =========================================================================
%
%  OUTPUTS
%    [m]            [mean PETH]
%    [t]               [times of PETH]
%
%
%  SEE ALSO
%
%   [PETH, DeconvolvePETH]
%
% [Ralitsa Todorova, Michaël Zugaro] [2018-2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%--------------------------------------------------------------------------

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



