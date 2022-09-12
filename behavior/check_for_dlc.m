function found_one = check_for_dlc(basepath,basename)
%check_for_dlc - folder checker for dlc files
%
%  Quick checker for Deep Lab Cut files within the provided basepath.
%
%
%  INPUTS
%    basepath       Full relevant folder destination as character string.
%    basename       Relevant basename to check as character string
% 
%  OUTPUTS
%    found_one      Array of logical values denoting whether or not DLC is
%                   correctly run for each session.
%
% AYA Lab 2022
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load(fullfile(basepath,[basename,'.MergePoints.events.mat']))

for k = 1:length(MergePoints.foldernames)
    dlc_flag(k) = isempty(dir(fullfile(basepath,MergePoints.foldernames{k},'*DLC*.csv')));
end
files = dir(basepath);
files = files(~contains({files.name},'Kilosort'),:);
dlc_flag(k+1) = isempty(dir(fullfile(files(1).folder,'*DLC*.csv')));

found_one = ~all(dlc_flag);
end