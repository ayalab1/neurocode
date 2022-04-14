function grouped = Group(varargin)

% Provide as many vectors as you want. They will be grouped in a matrix
% in the same format as spikes are grouped ([values id])
% vector1 1
% vector2 2
% ... and so on. Useful for creating grouped events.
% ATTENTION: rows are automatically sorted.
%
% Copyright (C) 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

grouped = [];
% Are we dealing with vectors?
if any(cellfun(@(x) min(size(x)), varargin)==1),
	vectors = true; else vectors = false; 
end
for i=1:length(varargin),
	if vectors
    grouped = [grouped; varargin{i}(:) i*ones(size(varargin{i}(:),1),1)];
	else
		try
		grouped = [grouped; varargin{i} i*ones(size(varargin{i}(:,1),1),1)];
		catch
			error('If you provide matrices, please make sure these matrices have the same number of columns.');
		end
	end
end

grouped = sortrows(grouped);