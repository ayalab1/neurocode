function filtered = FilterLFP(lfp,varargin)

%FilterLFP - Filter the local field potentials, e.g. in the theta band.
%
% Filter the local field potentials using a cheby2 filter.
%
%  USAGE
%
%    filtered = FilterLFP(lfp,<options>)
%
%    lfp            local field potentials <a href="matlab:help samples">samples</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'passband'    frequency range (default = [4 10], 'theta')
%     'order'       filter order (default = 5, but see NOTE below)
%     'ripple'      filter ripple (default = 20)
%     'filter'      choose filter type between 'cheby2' (default) and 'fir1'
%    =========================================================================
%
%  NOTE
%
%    The passband can be supplied either explicitly, e.g. [30 80], or by name,
%    by choosing among the following predefined frequency bands:
%
%        delta      [0 6] (order = 9)
%        theta      [4 10]
%        spindles   [9 17]
%        gamma      [30 80]
%        ripples    [100 250]

% Copyright (C) 2004-2021 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
passband = [4 10];
order = 5;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(lfp),
	error('Parameter ''lfp'' is not a matrix (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).');
end

% Parse parameter list
i = 1;
while i <= length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'passband',
			passband = varargin{i+1};
			if isastring(passband),
				switch(lower(passband)),
					case 'delta',
						passband = [0 6];
						order = 9;
					case 'theta',
						passband = [4 10];
					case 'spindles',
						passband = [9 17];
					case 'gamma',
						passband = [30 80];
					case 'ripples',
						passband = [100 250];
					otherwise,
						error(['Unknown frequency band ''' passband ''' (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).']);
				end
			elseif ~isdvector(passband,'#2','>=0'),
				error('Incorrect value for ''passband'' (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).');
			end
			varargin = {varargin{[1:(i-1) (i+2):length(varargin)]}};
		otherwise,
			i = i+2;
	end
end

%% Insertion by Raly
% original lfp
t0 = lfp(:,1);
ok = diff([-Inf;t0])>0;

% interpolate
t = (lfp(1,1):(1/1250):lfp(end,1))'; % interpolate the LFP to 0.8 ms bins, because the filter parameters work well for this sampling rate
interpolated = interp1(lfp(ok,1),lfp(ok,2),t);

if passband(1)>0
    slow = FilterLFP([t interpolated],'passband',[0 passband(1)],'order',order,varargin{:});
    if passband(2)==Inf
        filtered = [t interpolated-slow(:,2)];
    else
        filtered = FilterLFP([t interpolated-slow(:,2)],'passband',[0 passband(2)],'order',order,varargin{:});
    end
else
    filtered = Filter([t interpolated],'passband',passband,'order',order,varargin{:});
end

% interpolate back
filtered = [t0 interp1(filtered(:,1),filtered(:,2),t0)];
if isnan(filtered(end,2)), filtered(end,2) = filtered(end-1,2); end % make sure we don't lose the last bin due to interpolating above
filtered(any(isnan(filtered),2),:) = []; % remove any errors due to unforeseen nans

%% Before insertion it was just this line:
% filtered = Filter([t interpolated],'passband',passband,'order',order,varargin{:});


