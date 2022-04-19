function logical = IntervalsIntersect(intervals1, intervals2)

% This function returns a logical vector the size of intervals1, with ones 
% if any portion of intervals2 is present in the respective interval1.
% 
% This can be useful if, for example, you want to have the identity of
% the windows intervals1 which contain ripples (intervals2), without
% limiting this to the windows containing the ripple peak, but any part 
% of the ripple.
%
% Copyright (C) 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if isempty(intervals1) || isempty(intervals2), logical = false(size(intervals1,1),1); return; end

% start or end of intervals2 is included within an interval1. This should be counted.
l1 = CountInIntervals(intervals2(:,1), intervals1)>0;
l2 = CountInIntervals(intervals2(:,2), intervals1)>0;
% intervals2 completely incorporates ones of intervals1.
l3 = InIntervals(intervals1(:,2), intervals2);

logical = l1 | l2 | l3;
