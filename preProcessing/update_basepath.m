function update_basepath
%
% 
%   USAGE: update_basepath
%
%   After copying the files generated during preprocessing a given
%   session outside in your local drive (or somewhere else), run 
%   update_basepath inside the folder where you want your files to be and 
%   it will update the basename and basepath to match the current location
%
%   a copy of the original path where preprocessing was run is kept as 
%   'basepath_original' inside preprocessSession_params.mat
%
%
% [ayalab - 2025]

%%

%get final name to add to session files
basename = basenameFromBasepath(pwd);

%list all files to search through
fileList = dir(fullfile(pwd, '*.*'));

% Find session and session params to replace basepath
for k = 1:length(fileList)
    fileName = fileList(k).name;
    
    % Find and rewrite path for session.mat
    if contains(fileName, '.session.mat')
        
        tmp_basename = extractBefore(fileName,".session.mat");
        
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

% Find files that need correction and replace them
fileList = dir([tmp_basename '.*']);

for k = 1:length(fileList)
    filename_old = fileList(k).name;
    filename_new = strrep(filename_old, tmp_basename, basename);
    if ~exist(fullfile(pwd, filename_new),'file')
        movefile(filename_old, filename_new); 
    end
end

disp(['file names updated']);

end