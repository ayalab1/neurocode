%mmax - Return overall min and max values for any number of arrays.
% renamed mmax because Matlab2016 contains a different built-in function minmax
%
%  USAGE
%
%    [m,i] = minmax(x,y,...)
%
%    x,y,...        variables
%
%  OUTPUT
%
%    m              minima and maxima of x,y,... (one pair per line)
%    i              indices of minimum and maximum values

% Copyright (C) 2015 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [m,i] = mmax(varargin)

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help minmax">minmax</a>'' for details).');
end

for k = 1:length(varargin),
	v = varargin{k};
	[m(k,1),i(k,1)] = min(v(:));
	[m(k,2),i(k,2)] = max(v(:));
end
