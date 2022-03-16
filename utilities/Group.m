function grouped = Group(varargin)

% Provide as many vectors as you want. They will be grouped in a matrix
% in the same format as spikes are grouped:
% vector1 1
% vector2 2
% ... and so on. Useful for creating grouped events.
% ATTENTION: rows are automatically sorted.


%% Loop equivalent
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

% %% Vectorized version code (not actually faster)
% values = cell2mat(cellfun(@(x)cat(1,x(:)),varargin','UniformOutput',false));
% groups = repelem(1:length(varargin),cellfun(@numel,varargin));
% grouped = [values groups];
%
grouped = sortrows(grouped);