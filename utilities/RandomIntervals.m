function r = RandomIntervals(n, sumDuration)

%RandomIntervals - generates n random intervals from 0 to 1, covering a total of 'sumDuration'
%
% EXAMPLE
%
% Generate random surrogate ripple events 
% nRipples = size(ripples.timestamps,1);
% totalRippleDuration = sum(diff(ripples.timestamps,[],2));
% totalRecordingDuration = MergePoints.timestamps(end);
% p = totalRippleDuration/totalRecordingDuration; % the proportion of total time spent in ripples
% r = RandomIntervals(size(ripples,1),p);
% randomIntervals = r*totalRecordingDuration;
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~InIntervals(sumDuration,[0 1]),
    error('Overall duration of the intervals should be between 0 and 1');
end
r = zeros(n,2);
r(:,1) = sortrows(rand(n,1));
r(:,2) = r(:,1)+sumDuration.*([r(2:end,1); 1]-r(:,1)); % on avarage, that's "sumDuration"
r(1) = r(1)-sumDuration*(r(1));