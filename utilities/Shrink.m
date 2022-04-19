function shrunk = Shrink(matrix,rowBinSize,columnBinSize)

%Shrink - reduce matrix resolution by shrinking it (akin to pixelating)
% This returns a matrix 'shrunk' which is produced by taking the
% mean of [columnBinSize rowBinSize]-sized blocks within the input matrix.
% From this follows that columnBinSize needs to be a factor of size(matrix,2),
% and rowBinSize needs to be a factor of size(matrix,1). If not, nans will be
% be appended (as equally spaced as possible).
%
% Copyright (C) 2019 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


n0 = size(matrix);
n = ceil(n0./[rowBinSize columnBinSize]).*[rowBinSize columnBinSize];

% This is a heuristic to make the code run faster: 
if columnBinSize>rowBinSize, 
    shrunk = Shrink(matrix',columnBinSize,rowBinSize)';
    return;
end

if n0(1)~=n(1),
    if n0(1)==1, matrix = [matrix; nan(n(1)-1,n0(2))]; else
    nanRows = round(linspace(1,n(1),2+n(1)-n0(1))); nanRows([1 end]) = []; % Equally spaced rows of nans
    matrix0 = matrix;
    matrix = nan(n(1),n0(2));
    matrix(~ismember(1:n(1),nanRows),:) = matrix0;
    end
end

if n0(2)~=n(2),
    if n0(2)==1, matrix = [matrix nan(n(1),n(2)-1)]; else
	nanColumns = round(linspace(1,n(2),2+n(2)-n0(2))); nanColumns([1 end]) = []; % Equally spaced columns of nans
    matrix0 = matrix;
    matrix = nan(n);
    matrix(:,~ismember(1:n(2),nanColumns),:) = matrix0; % fill non-nan values
   end
end


if columnBinSize == 1 && rowBinSize == 1,
    shrunk = matrix;
    return
end

order = [];
for i=1:columnBinSize,
    if i<columnBinSize,
        order = [order find(rem(1:size(matrix,2),columnBinSize)==i)];
    else
        order = [order find(rem(1:size(matrix,2),columnBinSize)==0)];
    end
end

% This code was produced by trial and error. Too complicated to understand, let alone annotate.
shrunk = matrix(:,order);
shrunk = reshape(shrunk,rowBinSize,[]);
shrunk = reshape(shrunk,[],columnBinSize)';
shrunk = reshape(shrunk,1,[])';
shrunk = reshape(shrunk,columnBinSize*rowBinSize,[]);
shrunk = nanmean(shrunk);
shrunk = shrunk(:);
shrunk = reshape(shrunk,size(matrix)./[rowBinSize columnBinSize]);

