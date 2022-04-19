function s = semedian(x,dim)

%semedian - Compute standard error of the median.
%
%  USAGE
%
%    s = semedian(x,dim)
%
%    x              vector or matrix over which the error should be computed
%    dim            optionally, dimension along which to perform the analysis
%
% Copyright (C) 2013-2022 by Michaël Zugaro, Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help semedian">semedian</a>'' for details).');
end
if ~isdmatrix(x) && ~isdvector(x),
  error('Incorrect input - use vector or matrix (type ''help <a href="matlab:help semedian">semedian</a>'' for details).');
end

if ~exist('dim','var'), dim = 1; end
if ~isdvector(dim,'#1')
    error('Incorrect value for parameter ''dim'' (type ''help <a href="matlab:help semedian">semedian</a>'' for details).');
end

if any(size(x)==1), x = x(:); end

if exist('dim','var') && dim==2
    s = semedian(x')';
    return
end

n = size(x,1);
m = repmat(nanmedian(x),n,1);
s = sqrt( nansum((x-m).^2) / (n*(n-1)) );
