function Logical = Unfind(indices, n)
 
%
% When you have a logical vector, 'find' is useful as it gives you the
% indices of the non-zero values.
% This function performs the opposite operation, giving you a logical
% vector with ones in the positions of the indices provided.
%
% If you want the length of the resulting logical to be different from the
% last index provided, provide n (default = last index).
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
 
if nargin<2,
    n=max(indices);
end
 
Logical = false(n,1);
Logical(indices) = true;