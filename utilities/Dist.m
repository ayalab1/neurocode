function varargout = Dist(resolution, varargin)

% Dist - Outputs the normalised distributions of the data provided at a given resolution.
% [d,x] = Dist(resolution, data),
% where 'd' is the normalised distribution of the data for points 'x'.
% Dist(resolution, data) without output arguments plots said distributions.
%
% This is very similar to "hist" but it normalizes the distribitions and 
% accepts a variety of ways to input the data. 
%
% USAGE
%
% FOR VECTORS:
% [d,x] = Dist(100,vector1,vector2,...,vectorN);
% produces the pdf-s of the vectors provided, fitting in the same x-axis.
% The vectors provided need not have the same size.
% FOR A MATRIX:
% [d,x] = Dist(100,matrix) is equivalent to [d,x] = Dist(100,column1,column2,...,columnN);
% FOR GROUPED DATA:
% [d,x] = Dist(100,GroupedData,'grouped')
% will treat the last column of the matrix GroupedData as a grouping variable (from 1 to N),
% so it is equivalent to [d,x] = Dist(100,GroupedData(GroupedData(:,end)==1),GroupedData(GroupedData(:,end)==2),...)
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


barplot = 'off';
if strcmp(varargin{end}, 'bar')
    barplot = 'on';
    varargin = varargin(1:end-1);
end


if length(varargin)==1 && ~isvector(varargin{1})
    for i=2:size(varargin{1},2)
        varargin{i} = varargin{1}(:,i);
    end
    varargin{1} = varargin{1}(:,1);
elseif length(varargin)>1 && (strcmp(varargin{2}, 'group') || strcmp(varargin{2}, 'groups') || strcmp(varargin{2}, 'grouped')),
    for i=max(varargin{1}(:,end)):-1:1 %flipped so 1 comes last and we don't override the data
        varargin{i} = varargin{1}(varargin{1}(:,end)==i,1);
    end
    if (strcmp(varargin{2}, 'group') || strcmp(varargin{2}, 'groups') || strcmp(varargin{2}, 'grouped')), varargin(2) = []; end
else % treat each input as a separate variable. Make sure they are vertical vectors:
    for i=1:length(varargin), varargin{i} = varargin{i}(:); end
end

n = length(varargin);

for i=1:n
    if isempty(varargin{i})
        minmaxima(i,1:2) = nan;
    else
        minmaxima(i,1) = min(varargin{i}(:));
        minmaxima(i,2) = max(varargin{i}(:));
    end
end

limits = [min(minmaxima(:)) max(minmaxima(:))];
limits = [limits(1)-diff(limits)/5 limits(2)+diff(limits)/5]; % broaden by 40%
[~, t] = hist(limits, resolution);

for i=1:n
    [h(:,i), ~] = hist(varargin{i}, t);
end

h = h./repmat(sum(h), size(h,1), 1);

if nargout==0
    
    if strcmp(barplot, 'on')
        handle = bar(t,h);
    else
        handle = plot(t, h);
    end
    
    legend(num2str([1:n]'));
end
if nargout>0
    varargout{1} = h;
    varargout{2} = t;
end