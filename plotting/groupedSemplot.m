function groupedSemplot(g,y,varargin)

%groupedSemplot - alternative call to <a href="matlab:help semplot">semplot</a>
% You want to call semplot but your data is not in a matrix form;
% Rather, you have a different number of observations for each 
% x-value. Provide the groups (g= binned x-values) and the function
% will call semplot for you.
%
% Copyright (C) 2020-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

bad = isnan(g); g(bad) = []; y(bad) = [];
[u,~,uu] = unique(g(:));

n = sum(g==mode(g));
matrix = nan(length(u),n);
rankG = zeros(size(g));

for i=u'
	rankG(g==i) = cumsum(g(g==i)==i);
end

matrix(sub2ind(size(matrix),uu,rankG)) = y;

x = unique(g);
semplot(x,matrix,varargin{:});