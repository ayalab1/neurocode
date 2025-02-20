function files = recursive_dir(rootDir, pattern, excludeDirs)
%RECURSIVE_DIR Recursively search for files matching a pattern, excluding specified directories.
%   FILES = RECURSIVE_DIR(ROOTDIR, PATTERN, EXCLUDEDIRS) searches for files
%   matching PATTERN within ROOTDIR and all its subdirectories, excluding
%   directories listed in EXCLUDEDIRS.
%
%   Inputs:
%       ROOTDIR    - A string specifying the root directory to start the search.
%       PATTERN    - A string specifying the file pattern to match (e.g., '*.avi').
%       EXCLUDEDIRS - A cell array of strings specifying folder names to exclude
%                    from the search (e.g., {'.phy', 'temp'}).
%
%   Output:
%       FILES      - A cell array of strings containing the full paths to all
%                    matching files.
%
%   Example:
%       % Find all .avi files in 'U:\data\hpc_ctx_project', excluding '.phy' folders
%       files = recursive_dir('U:\data\hpc_ctx_project', '*.avi', {'.phy'});
%
%   See also DIR, FULLFILE.


% Initialize the list of files
files = [];

% Get the list of files and folders in the current directory
dirData = dir(fullfile(rootDir, pattern));

% Filter out directories from the list
dirData = dirData(~[dirData.isdir]);

% Add the full path to the files
files = [files; fullfile({dirData.folder}, {dirData.name})'];

% Get the list of subdirectories
subDirs = dir(rootDir);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name}, {'.', '..'}));

% Recursively call this function on each subdirectory, excluding specified folders
for i = 1:length(subDirs)
    subDirPath = fullfile(rootDir, subDirs(i).name);

    % Skip excluded directories
    if ~any(strcmp(subDirs(i).name, excludeDirs))
        files = [files; recursive_dir(subDirPath, pattern, excludeDirs)];
    end
end
end