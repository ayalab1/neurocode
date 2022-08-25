function [catted] = catOverlap(unCatted)

%catOverlap - concatenate overlapping timebins
%
%  Function for concatenating event periods that overlap
%
%  USAGE
%
%    % Combining event bins for differing events
%    [catted] = catOverlap(sort([ripples.timestamps; HSE.timestamps]));
%
%  INPUTS
%    unCatted       Nx2 matrix of start and stop times for N events,
%                   ideally sorted by start time
% 
%  OUTPUTS
%    catted         sorted Mx2 matrix, where M <= N (event
%                   times will be concatenated into one event if 
%                   overlapping). 
%
%  NOTE
%
%    Input times do not need to be sorted, but must follow the format of
%    column 1 being start times and column 2 being end times
%
%  EXAMPLES
%
%    [ie. “catOverlap([1 2; 1.5 7; 8 9]);  % returns [1 7; 8 9]”]
%    [ie. “catOverlap([1.5 7; 8 9; 1 2]);  % returns [1 7; 8 9]”]
%
%  SEE
%
%    

% Copyright (C) 2022 by Lindsay Karaba
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

[sorted, sI] = sort(unCatted(:,1));
sorted(:,2) = unCatted(sI,2);
catted = sorted(1,:);
cat_ct = 1;
for i = 2:size(sorted,1)
    if sorted(i,1) <= catted(cat_ct,2) %if start of next is less than end of last
        if catted(cat_ct,2) < sorted(i,2)
            catted(cat_ct,2) = sorted(i,2); %extend end evt time if needed
        end
    else
        cat_ct = cat_ct+1;
        catted(cat_ct,:) = sorted(i,:);
    end
end
end
