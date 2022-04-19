function y = linspaceVector(d1, d2)

%linspaceVector - produce linearly spaced vector with multiple [start stop] points. 
% It's like linspace, but d1 and d2 and vectors
% It's a faster equivalent to a loop calling linscape for each [d1 d2] pair:
% for i=1:length(d1),y = [y;linspace(d1(i),d2(i))']; end
%
% Example:
% linspaceVector([2;5;23],[2;8;23]) = [2;5;6;7;8;23]
%
% Copyright (C) 2019 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

n = d2-d1;

% repeat original value n times
y0 = repelem(d1,n+1);
nRows = size(y0,1);

% add 1 for each additional value after d1
indicesOfOriginalValues = cumsum([0;n(1:end-1)]+1);
originalValues = Unfind(indicesOfOriginalValues,nRows);
toAdd = CumSum(ones(nRows,1),originalValues)-1;

y = y0 + toAdd;
