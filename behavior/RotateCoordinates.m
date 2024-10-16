function [rotated] = RotateCoordinates(xy,angle,center)

% Rotates x and y coordinates around a center point by an angle.
% 
% (angle should be in radians)
% This function rotates x and y coordinates (e.g. for the positions file)
% around a center point (default=origin) by an angle.
%
%
%  USAGE
%
%    [xy] = RotateCoordinates(xy,angle,center)
%
%    xy             a matrix containing the [x y] coordinates to be rotated
%    angle          the angle of desired rotation (in radians)
%    center         the [x y] coordinates of the center of rotation
%
%  OUTPUT
%
%    rotated        modifed [x y] coordinates after the rotation
%
%
% EXAMPLE:
% positions(:,2:3) = RotateCoordinates(positins(:,2:3),pi/8,[50 60])
%
% Copyright (C) 2018-2024 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('center','var'),
    center = [0 0];
end

R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
rotated = bsxfun(@minus,xy,center);
rotated = rotated*R;
rotated = bsxfun(@plus,rotated,center);
