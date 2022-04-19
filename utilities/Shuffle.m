function shuffled = Shuffle(matrix,varargin)

%Shuffle - shuffle the elements of a matrix
% Independently shuffle the elements within each row (default) or column, or
% all elements ignoring the matrix dimensions. Optionally, you may select
% certain elements to remain fixed.
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'keep'        elements to keep in place (all other elements will be shuffled)
%     'dimension'   dimension along which to shuffle: 1 would shuffle columns,
%                   2, rows, ':' will shuffle all the values randomly
%                   (default = 2)
%    ===========================================================================
%
% OUTPUTS
%
%    shuffled       the shuffled version of 'matrix'
%
% Copyright (C) 2016-2022 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

keep = false(size(matrix));
dimension = 2;
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Shuffle">Shuffle</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'keep',
            keep = varargin{i+1};
            try 
                keep = logical(keep);
                if ~any(size(keep)==size(matrix)), error(''); end
            catch
                error('Incorrect value for property ''keep'' (type ''help <a href="matlab:help Shuffle">Shuffle</a>'' for details).');
            end
        case 'dimension',
            dimension = varargin{i+1};
            if (~isvector(dimension) || length(dimension) ~= 1) && ~isastring(dimension,':'),
                error('Incorrect value for property ''dimension'' (type ''help <a href="matlab:help Shuffle">Shuffle</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Shuffle">Shuffle</a>'' for details).']);
    end
end

if ~isdmatrix(matrix) && ~isdvector(matrix),
    error('Please provide a matrix of doubles to shuffle)');
end
if isempty(matrix)
    shuffled = matrix;
    return
end

if size(keep,1)~=size(matrix,1), keep = repmat(keep,size(matrix,1),1); end
if size(keep,2)~=size(matrix,2), keep = repmat(keep,1,size(matrix,1)); end

if dimension==1, shuffled = Shuffle(matrix','keep',keep','dimension',2)'; return; end
if isastring(dimension,':'), shuffled = matrix; shuffled(:) = Shuffle(matrix(:)','keep',keep(:)'); return; end

if any(keep(:)),
    [~,indices] = sort(rand(size(matrix')).*(~keep'));
else
    [~,indices] = sort(rand(size(matrix')));
end

shuffled = matrix(sub2ind(size(matrix),meshgrid(1:size(matrix,1),1:size(matrix,2))',indices'));
