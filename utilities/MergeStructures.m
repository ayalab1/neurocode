function merged = MergeStructures(structure1,structure2,dim)

% This only works if the two structures have the same fields and all the
% sizes are compatible

if ~exist('dim','var')
    dim = 1;
end

names = fieldnames(structure2);
for i=1:length(names)            
    if ~isstruct(structure1.(names{i}))
        merged.(names{i})= cat(dim,structure1.(names{i}), structure2.(names{i}));
    else % allow nesting
        merged.(names{i})= MergeStructures(structure1.(names{i}), structure2.(names{i}),dim);
    end
end

