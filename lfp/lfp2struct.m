function structure = lfp2struct(lfp,lfpInfo)
%lfp2struct: temporary function transforming lfp into structure format
%
%This function takes the outputs of getLFP (lfp data in matrix format and 
%additional metadata in lfpInfo) and stores them into a structure. This can
%be useful when calling not-yet-updated neurocode functions which require
%the structure format.
%
%%OUTPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels
%
% EXAMPLE
% [lfp,info] = getLFP(1);
% CFCPhaseAmp(lfp2struct(lfp,info)); % because CFCPhaseAmp is not yet updated and still requires structure input
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

structure = lfpInfo;
structure.timestamps = lfp(:,1);
structure.data = lfp(:,2:end);
