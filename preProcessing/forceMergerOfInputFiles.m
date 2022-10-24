function forceMergerOfInputFiles(basepath)

% [This function will produce a digitalIn.dat and analogin.dat files even if
%  they are missing from some of the subsessions]
%
%
%  INPUTS
%    [basepath]      [basepath of input files to be merged]
%
%
%  OUTPUTS
%    NA
%
%
% SEE ALSO
%
% [Ralitsa Todorova] [2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

cd(basepath);
basename = basenameFromBasepath(basepath);

if ~exist('digitalin.dat','file')
    list = dir([basename '*']); list(~cat(1,list.isdir),:) = [];
    isdig = false(size(list,1),1);
    for i=1:length(list)
        isdig(i) = exist(fullfile(list(i).name,'digitalin.dat'),'file');
    end
    if any(isdig)
        % Make a digitalin.dat file
        f = fopen('digitalin.dat','w');
        for i=1:length(list)
            if isdig(i)
                file = memmapfile(fullfile(list(i).name,'digitalin.dat'),'Format','int16','Writable',false);
                fwrite(f,file.data(:),'int16');
            else
                % make an empty binary of the correct length as what the binary would have been, which should be
                % half the size of the 'time.dat' file
                file = memmapfile(fullfile(list(i).name,'time.dat'),'Format','int16','Writable',false);
                data = reshape(file.data,2,[]); data(:) = 0;
                fwrite(f,data(1,:)','int16');
            end
        end
    end
end
fclose('all');

if ~exist('analogin.dat','file')
    list = dir([basename '*']); list(~cat(1,list.isdir),:) = [];
    isdig = false(size(list,1),1);
    for i=1:length(list)
        isdig(i) = exist(fullfile(list(i).name,'analogin.dat'),'file');
    end
    if any(isdig)
        % Make a digitalin.dat file
        f = fopen('analogin.dat','w');
        for i=1:length(list)
            if isdig(i)
                file = memmapfile(fullfile(list(i).name,'analogin.dat'),'Format','int16','Writable',false);
                fwrite(f,file.data(:),'int16');
            else
                % make an empty binary of the correct length as what the binary would have been, which should be
                % half the size of the 'time.dat' file
                file = memmapfile(fullfile(list(i).name,'time.dat'),'Format','int16','Writable',false);
                data = reshape(file.data,2,[]); data(:) = 0;
                fwrite(f,data(1,:)','int16');
            end
        end
    end
end
fclose('all');
