function EquateScales(varargin)

%EquateScales - Equalize the the x and y scales on all the axes of the current figure
% EquateScales('x') equates the x-axes only, and EquateScales('y') - the y-axes only.
%
% Copyright (C) 2020 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

axes = findall(gcf, 'type', 'axes');

normy = true;
normx = true;

if ~isempty(varargin),
    theseaxes = [];
    for i=1:length(varargin),
        ax = varargin{i};
        if strcmp(ax,'x'),normy=false;
        elseif strcmp(ax,'y'), normx=false;
        elseif isvector(ax),
            theseaxes = [theseaxes numel(axes)+1 - ax];
        end
    end        
else
    theseaxes = 1:numel(axes);
end

y = nan(numel(axes),2);
x = nan(numel(axes),2);
for i=theseaxes,
    axis = axes(i);
    y(i,:) = get(axis, 'YLim');
    x(i,:) = get(axis, 'XLim');
end

y = [min(y(:,1)) max(y(:,2))];
x = [min(x(:,1)) max(x(:,2))];

for i=theseaxes,
    axis = axes(i);
	if normy,set(axis, 'YLim', y);end
    if normx,set(axis, 'XLim', x);end
end
