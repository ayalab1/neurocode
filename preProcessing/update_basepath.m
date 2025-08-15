function update_basepath
%
%  [update_basepath] - [update basepath to match current path]
% 
%   Update basepath inside session files to match current path
%   (to correct for when data is preprocessed in a different drive)
%
%   session.general.basepath will be updated to new path
%   preprocessSession_params.mat will be updated to new path but will keep
%   a copy of original path as 'basepath_original'
%
%
% [ayalab - 2025]

%%

%get final name to add to session files
basename = basenameFromBasepath(pwd);

%list all files to search through
fileList = dir(fullfile(pwd, '*.*'));

% Find files that need correction and replace basepath
for k = 1:length(fileList)
    fileName = fileList(k).name;
    
    % Find and rewrite path for session.mat
    if contains(fileName, '.session.mat')
        fullFilePath = fullfile(pwd, fileName);
        
        % Load variables from the file
        load(fullFilePath);
        
        % Overwrite basepath
        session.general.basePath=pwd;
        session.general.name=basename;
        save(fullfile(pwd, [basename, '.session.mat']), 'session');
        
        %getting rid of the old one in case it had a different name
        if ~strcmp(fullfile(pwd, [basename, '.session.mat']),fullFilePath)
            delete(fullFilePath);
        end
    end
    
    % Find and rewrite path for preprocessSession_params.mat
    if contains(fileName, 'preprocessSession_params.mat')
        fullFilePath = fullfile(pwd, fileName);
        
        % Load variables
        load(fullFilePath);
        
        % Add info about new path (keeping old path for whatever)
        results.basepath_original=results.basepath;
        results.basepath=pwd;
        save(fullfile(pwd, ['preprocessSession_params.mat']), 'results');
    end
end

disp(['basepath updated']);

end