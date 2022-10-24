function lfp = GetAyaLFP(channels,varargin)

%[GetAyaLFP] - temporary, now obsolete, wrapper which calls getLFP
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

lfp = getLFP(channels,varargin{:});
% lfpstruct = getLFP(channels+1,varargin{:});
% try lfp = [lfpstruct.timestamps]; 
%     lfp(:,(1:size(lfpstruct.data,2))+1) = lfpstruct.data;
% catch
%     keyboard
% end
