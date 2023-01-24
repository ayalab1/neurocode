function merged = MergeStructures(structure1,structure2)

% This only works if the two structures have the same fields and all the
% sizes are compatible

names = fieldnames(structure2);
for i=1:length(names)              
    merged.(names{i})= cat(1,structure1.(names{i}), structure2.(names{i}));
end

