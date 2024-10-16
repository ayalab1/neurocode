function find_and_move_videos(basepath,local_drive_path)
% Matlab wrapper for find_and_move_videos.py
%
% find_and_move_videos will copy videos to dlc video path and also copy
% results back to original folder locations
%
% Input: 
%       basepath: session or project path
%       local_drive_path: local path where you dlc videos were live
%
% Example:
% 
%   basepath = 'Z:\Data\Can\OLM21'
%   local_drive_path = 'C:\Users\Cornell\dlc_videos_v2'
%   find_and_move_videos(basepath,local_drive_path)
%
% if python not detected, give your python path like so:
%                   run pyversion('C:\Users\Cornell\anaconda3\python.exe')
% 
% Ryan Harvey 2022

disp(pyversion)
pathToCode = fileparts(which('find_and_move_videos.py'));
if count(py.sys.path,pathToCode)==0
    % adds the code path to the python path
    insert(py.sys.path,int32(0),pathToCode) 
end
pyOut = py.find_and_move_videos.main(basepath,local_drive_path);
end