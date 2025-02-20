function varargout = semplot(x,y,color,smooth,solid)

%semplot - plot mean (line) +/- s.e.m. (shaded area) of a matrix "y"
% semplot(x,y,color,smooth)
%
% Copyright (C) 2017 by Ralitsa Todorova
%
% SEE ALSO shadedErrorBar 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('solid','var')
    solid = false;
end

if nargin<2 || (nargin==2 && ischar(y))
    if nargin==2
        color = y;
    end
    y = x;
    x = 1:size(y,2);
end

if ~exist('color','var')
    color = [0 0 0];
end

if ~exist('smooth','var')
    smooth = 0;
end

if size(y,2)~=length(x)
    y = y';
    if size(y,2)~=length(x)
        error('Y should have one column for each element in X');
    end
end

if isvector(y)
    handles = plot(x,Smooth(y,smooth),'color',color);
    if nargout>0, varargout{1} = handles; end
    hold on;
    return
end

bad = isnan(nanmean(y));
if any(bad)
    semplot(x(~bad),y(:,~bad),color,smooth);
end

xx = [x(:);flipud(x(:))];
yy = [Smooth(nanmean(y)'-nansem(y)',smooth); Smooth(flipud(nanmean(y)'+nansem(y)'),smooth)];
y = Smooth(nanmean(y),smooth);
handles = fill(xx,yy,color);

if solid
    set(handles,'FaceColor',mean([color;1 1 1]),'edgeAlpha',0);
else % transparent
    set(handles,'FaceAlpha',0.5,'edgeAlpha',0);
end

hold on;
plot(x,y,'color',color,'linewidth',2);

if nargout>0, varargout{1} = handles; end