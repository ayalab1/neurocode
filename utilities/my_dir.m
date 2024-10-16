function fileList = my_dir(directory)
% my_dir Returns a table of file and directory information for a specified directory.
%
%   fileList = MY_DIR(directory) returns a table containing information about
%   the files and directories in the specified directory. The input argument
%   'directory' is a string specifying the directory to be listed. The output
%   'fileList' is a table where each row represents a file or directory in the
%   specified directory, and the columns represent the different properties of
%   each file or directory.
%   
%   I made this because the default disp output is not visable. Also, the
%   default struct output is not fun to work with. 
%
%   Example:
%       fileList = my_dir('/path/to/your/directory');
%
%   See also: DIR, STRUCT2TABLE.
%   
%   Ryan H
%
% Get the directory listing using dir
dirStruct = dir(directory);

% Convert the directory structure to a table
fileList = struct2table(dirStruct);
end
