function [strength,dN] = ReactivationStrength(spikes,templates,varargin)

%ReactivationStrength - Assess reactivation strength of cell assemblies.
%
% Estimates how similar spiking activity is to given templates, across time.
% Time bins can be automatically determined using a fixed bin size, or provided
% as an explicit list (e.g. computed using theta phases).
%
%  USAGE
%
%    strength = ReactivationStrength(spikes,templates,<options>)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bins'        list of [start stop] for all bins
%     'binSize'     bin size in s (default = 0.050)
%     'overlap'     overlap between successive bins (default = binSize/2)
%     'step'        step between successive bins (default = binSize/2)
%    =========================================================================
%
%  OUTPUT
%
%    strength       reactivation strength across time (if time bins are not
%                   explicitly provided, the first bin is centered on the
%                   first spike)
%
%  SEE
%
%    See also ActivityTemplates.

% Copyright (C) 2016-2018 by Michaël Zugaro, Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
defaultBinSize = 0.050;
binSize = [];
overlap = [];
step = [];
bins = [];

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'binsize',
            binSize = varargin{i+1};
            if ~isdscalar(binSize,'>0'),
                error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).');
            end
        case 'overlap',
            overlap = varargin{i+1};
            if ~isdscalar(overlap,'>0'),
                error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).');
            end
        case 'step',
            step = varargin{i+1};
            if ~isdscalar(step,'>0'),
                error('Incorrect value for property ''step'' (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).');
            end
        case 'bins',
            bins = varargin{i+1};
            if ~isdmatrix(bins,'@2'),
                error('Incorrect value for property ''bins'' (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).']);
    end
end

% Option bins is incompatible with binSize, step and overlap
if ~isempty(bins) && (~isempty(binSize) || ~isempty(overlap) || ~isempty(step)) ,
	error('Parameter ''bins'' is incompatible with ''binSize'', ''overlap'' and ''step'' (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).');
end
% Unless explicit bins were provided, update binSize and step/overlap
if isempty(bins),
	if isempty(binSize),
		binSize = defaultBinSize;
	end
	if isempty(step) && isempty(overlap),
		step = binSize/2;
	else
		if isempty(step),
			step = binSize-overlap;
		elseif step ~= binSize-overlap,
			error('Incompatible ''step'' and ''overlap'' parameters (type ''help <a href="matlab:help ReactivationStrength">ReactivationStrength</a>'' for details).');
		end
	end
end

if isempty(templates)
    strength = zeros(0,1); dN = zeros(0,1); return
end

% Shift spikes to start at 0 and bin them
nUnits = size(templates,2);
spikes = spikes(spikes(:,2)<=nUnits,:);
spikes = sortrows(spikes,1);
id = spikes(:,2);
if isempty(bins),
	spikes(:,1) = spikes(:,1);
	bins = (spikes(1,1):step:(spikes(end,1)-binSize))';
	bins(:,2) = bins + binSize;
end

% Create and compute spike count matrix
nBins = size(bins,1);
dN = zeros(nBins,nUnits);
for unit = 1:nUnits,
	dN(:,unit) = CountInIntervals(spikes(id==unit,1),bins);
end
%  dN = squeeze(reshape(dN,1,[],nUnits));
dN = zscore(dN);

% Compute reactivation strengths
nTemplates = size(templates,3);
strength = zeros(nBins,nTemplates);
for i = 1:nTemplates,
    template = templates(:,:,i);
    % Set the diagonal to zero to not count coactivation of i and j when i=j
    template = template - diag(diag(template));
    strength(:,i) = nansum(dN*(template).*dN,2);
end
t = nanmean(bins,2);
strength = [t strength];

