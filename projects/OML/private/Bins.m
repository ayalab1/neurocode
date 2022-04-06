function bins = Bins(startpoint,endpoint,window,step)

% Just a quick function to make bins spanning from
% t=startpoint to t=endpoint, with a given window and step
% Copyright (C) 2019 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin<4,
	if nargin<3,
		window = 1;
	end
	step = window;
end

bins = startpoint:step:endpoint-window;
bins = bins(:);
bins = [bins bins+window];