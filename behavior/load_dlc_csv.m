function df = load_dlc_csv(file_path)
% load_dlc_csv: general function to load dlc csv files
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
% TODO: Currently, multi-tracking points have only basic column name 
%           (x,y,likelihood), but more detailed info will be needed such 
%           as bodypart, animal, etc.
%
% Ryan H 2022

header_i = 1;
while true
    opts = detectImportOptions(file_path,'NumHeaderLines',header_i);
    df = readtable(file_path,opts);
    % get names of fields, these will be as long as tracker points
    % used times 3 because [x,y,likelihood]
    field_names = df.Properties.VariableNames;
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

% make sure all columns are doubles (known issue where some were strings)
for fn = field_names
    if ~isa(df.(fn{1}),'double')
        df.(fn{1}) = str2double(df.(fn{1}));
    end
end
end