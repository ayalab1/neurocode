function interval = FindInterval(logical)

% Copyright (C) 2016 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

logical = [logical(:)' 0];
starts = strfind(logical,[0 1])+1;
stops = strfind(logical,[1 0]);
if isempty(starts),
    if isempty(stops),
        if logical(1),
            interval = [1 length(logical)];
        else
            interval = [];
        end
        return
    else
        starts = 1;
    end
end

if isempty(stops) || starts(end)>stops(end),
    stops = [stops length(logical)];
end

if starts(1)>stops(1),
    starts = [1 starts];
end

interval = [starts' stops'];