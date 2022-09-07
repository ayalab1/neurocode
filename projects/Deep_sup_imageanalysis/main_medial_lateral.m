%% Go over all files in current folder and compute cell depths (output will be mat files! existing files will be rewritten)
close all

% inpath = "D:\微云同步盘\George's Cloud\Abroad\NYU\工作\Buzsaki Lab\Projects\Roman\Development\Histology analysis\ImageJ results\";
files = dir(fullfile(strcat(pwd,'\*.mat')));
if isempty(files)
    error("please run main first to get cell depth and everything!")
end
nFiles = length(files);
for i=1:nFiles
    [filepath,name,ext] = fileparts(files(i).name);
    name
    [medial_lateral_position, normalized_medial_lateral_position] = medial_lateral_distance(name);
end
% fclose('all')
