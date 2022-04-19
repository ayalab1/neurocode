function [xy] = RotateCoordinates(xy,angle,center)

% angle should be in radians
% This function rotates x and y coordinates (e.g. for the positions file)
% around a center point (default=origin) by an angle.
% EXAMPLE:
% positions(:,2:3) = RotateCoordinates(positins(:,2:3),pi/8,[50 60])
%
% Copyright (C) 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('center','var'),
    center = [0 0];
end

R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
xy = bsxfun(@minus,xy,center);
xy = xy*R;
xy = bsxfun(@plus,xy,center);
