function merged = MergeStructures(structure1,structure2,dim)

%MergeStructures - merge two structures into a single structure 
% 
% This function concatenates each of the fields of structure1 and 
% structure2 along dimention 'dim' (default = 1).
% If the sizes of the field contents are not compatible (e.g. cannot 
% concatenate two matrices of different sizes along the other dimension), 
% then the contents of the fields of structure1 and structure2 are stacked
% in a cell in the merged structure.
%
% Copyright (C) 2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('dim','var')
    dim = 1;
end

names = fieldnames(structure1);
for i=1:length(names)
    % Determine that the same field exists in structure2:
    if ~isfield(structure2,names{i})
        warning(['The second structure does not contain a field named ''' names{i} ''' Copying the field from the first structure only...']);
        merged.(names{i}) = structure1.(names{i});
        continue
    end
    % Determine the type of field we are dealing with:
    fieldType = 'concatenatable';
    if isstruct(structure1.(names{i})); fieldType = 'struct'; end
    if iscell(structure1.(names{i})); fieldType = 'cell'; end
    switch fieldType
        case 'cell'
            merged.(names{i}) = helperCellFields(structure1,structure2,names{i},dim);
        case 'concatenatable'
            % general case: attempt to concatenate the two fields:
            try
                merged.(names{i}) = cat(dim,structure1.(names{i}), structure2.(names{i}));
            catch
                warning(['The sizes of the ''' names{i} ''' field in the two structures don''t fit. Converting into cell...']);
                merged.(names{i}) = helperCellFields(structure1,structure2,names{i},dim);
            end
        case 'struct' % allow nesting: merge the fields within the fields
            merged.(names{i}) = MergeStructures(structure1.(names{i}), structure2.(names{i}),dim);
    end
end

function cellOutput = helperCellFields(structure1,structure2,name,dim)

if ~iscell(structure1.(name))
    structure1.(name) = {structure1.(name)};
end
if ~iscell(structure2.(name))
    structure2.(name) = {structure2.(name)};
end
cellOutput = cat(dim,structure1.(name), structure2.(name));