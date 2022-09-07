%% Go over all files in current folder and compute cell depths (output will be mat files! existing files will be rewritten)
close all
disp('Remember to check if the scale is correct!')
% inpath = "D:\微云同步盘\George's Cloud\Abroad\NYU\工作\Buzsaki Lab\Projects\Roman\Development\Histology analysis\ImageJ results\";
files = dir(fullfile(strcat(pwd,'\*.csv')));
nFiles = length(files);
for i=1:nFiles
    [filepath,name,ext] = fileparts(files(i).name);
    name
    [pFoot,sign,depth_abs,data_X,data_Y,t_a,line_X_down,line_Y_down,line_X_up,line_Y_up,file,scale,SP_thickness,depth_norm] = cell_depth(name);
end
fclose('all')

warning('please manually save the cell depth file!')