function SaveCustomEvents(filename,events,names,varargin)

%SaveCustomEvents - Save events (ripple, deltas etc)
%
%  USAGE
%
%    SaveCustomEvents(filename,ripples,eventNames,options)
%
%    filename       file to save to
%    events         event timestamps
%    names          event description (1 for each column)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'overwrite'   overwrite file if it exists (default = 'on')
%    =========================================================================
%
%  SEE
%
%    See also SaveRippleEvents
%
% Copyright (C) 2019 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
overwrite = 'on';

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SaveCustomEvents">SaveCustomEvents</a>'' for details).');
end

for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SaveCustomEvents">SaveCustomEvents</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'overwrite',
			overwrite = varargin{i+1};
			if ~isstring(overwrite,'on','off'),
				error('Incorrect value for property ''overwrite'' (type ''help <a href="matlab:help SaveCustomEvents">SaveCustomEvents</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SaveCustomEvents">SaveCustomEvents</a>'' for details).']);
	end
end

if ~iscell(names),
    if size(events,2)==1, names = {names}; 
    else
        error('Incorrect value for input ''names''. Please provide a cell array input.');
    end
end

nCols = length(names);
n = size(events,1);
output.time = reshape(events(:,1:nCols)',[],1);
output.description = repmat(names(:),n,1);

if strcmp('overwrite','on'),
    SaveEvents(filename,output,'overwrite','on');
else
    SaveEvents(filename,output);
end
