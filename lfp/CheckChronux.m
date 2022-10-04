%CheckChronux - Make sure chronux is installed and functional.
%
%  USAGE
%
%    CheckChronux

% Copyright (C) 2007-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function CheckChronux(f)

if nargin == 0,
	f = 'chronux';
	message = 'This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.';
else
	message = ['This function requires the <a href="http://www.chronux.org">chronux</a> toolbox (could not find ''' f ''').'];
end

if isempty(which(f)),
	error(message);
end

