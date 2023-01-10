function overlap = OverlapIntervals(varargin)

% This is a helper function which will return the overlap of all provided arguments.
%
% EXAMPLE:
% overlap = OverlapIntervals(ripples,preSleep, sws); 
% would only return intervals contained in "ripples", "preSleep" as well as "sws"
%
% Copyright (C) 2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

overlap = [-Inf Inf]; % start with an interval encompassing everything
for i=1:length(varargin)
    % everything outside of the current intervals should be excluded from the overlap:
    outside = SubtractIntervals([-Inf Inf],varargin{i}); 
    overlap = SubtractIntervals(overlap,outside);
end