function clims(varargin)

%clims - Equalize the the color limits on all the axes of the current figure
% Copyright (C) 2017 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

axes = findall(gcf, 'type', 'axes');

cLim = [];
theseaxes = [];

for i=1:length(varargin),
    if length(varargin{i})==2,
        cLim = varargin{i};
    elseif length(varargin{i})==1,
        theseaxes = [theseaxes (numel(axes)+1-varargin{i})];
    end
end

if isempty(theseaxes),
    theseaxes = 1:numel(axes);
end

if isempty(cLim),
    c = nan(numel(axes),2);
    for i=theseaxes,
        axis = axes(i);
        c(i,:) = get(axis, 'cLim');
    end
    cLim = [min(c(:,1)) max(c(:,2))];
end

for i=theseaxes,
    axis = axes(i);
    set(axis, 'clim', cLim);
    set(axis, 'clim', cLim);
end
