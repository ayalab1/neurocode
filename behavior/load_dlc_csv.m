function df = load_dlc_csv(file_path)
% 
% General function to load dlc csv files
%
% works by iterating over headers until x,y,likelihood are found, then
% makes sure all columns are type double. 
% 
% Iterating over headers are needed as multi-animal tracking has more
% headers than the standard single animal tracking
%
% Input: 
%       file_path: fullpath+name of csv file
% Output:
%       df: table containing tracking results [x,y,likelihood] for each
%       tracking point. 
%
%
% Ryan H 2022

unique_fields = {};
header_i = 1;
while true
    opts = detectImportOptions(file_path,'NumHeaderLines',header_i);
    df = readtable(file_path,opts);
    % get names of fields, these will be as long as tracker points
    % used times 3 because [x,y,likelihood]
    field_names = cellfun(@(x) regexprep(x,'_\d',''),df.Properties.VariableNames,'UniformOutput',false);
    unique_fields(header_i,:) = field_names;
    if any(contains(field_names,'x')) &&...
            any(contains(field_names,'y')) &&...
            any(contains(field_names,'likelihood'))
        break
    end
    header_i = header_i + 1;
    
    % error if x,y,likelihood not found
    if header_i > 50
        error([file_path,' not a dlc file'])
    end
end

% rename fields with concat headers offset 1 to detected [x,y,likelihood]
fields_ = {};
for i = 1:size(unique_fields,2)  
    fields_{i} = strjoin(unique_fields(:,i),'_');
end
df.Properties.VariableNames = fields_;

% make sure all columns are doubles (known issue where some were strings)
for fn = fields_
    if ~isa(df.(fn{1}),'double')
        df.(fn{1}) = str2double(df.(fn{1}));
    end
end
end