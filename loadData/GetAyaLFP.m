function lfp = GetAyaLFP(channels,varargin)

% Loads the lfp saved in CellExplorer format into a structure-free matrix format, 
% which FMAToolbox users might find useful. The channel should be provided in
% 0-format (same as neuroscope). 
%
%  USAGE
%
%    [spikes,regionID,regionNames] = GetAyaLFP(channel,<options>);
%
%    channels       list of channels to load
%    <options>      optional list of property-value pairs that will be passed
%                   on to the ayalab "getLFP" (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%    'basepath'     folder in which .lfp file will be found (default = pwd)
%    'basename'     file name to load (in case of multiple sessions in same folder)
%    'intervals'    list of time intervals [0 10; 20 30] to read from 
%                   the LFP file (default is [0 inf])
%    'downsample'   factor to downsample the LFP by ('downsample' = 5 would load 
%                   a 1250Hz .lfp file at 250Hz; default = 1);
%    'noPrompts'    boolean to supress any user prompts (default = true)
%    'fromDat'      option to load directly from .dat file (default = false) instead
%                   of the .lfp file
%    =========================================================================
%
%  OUTPUT
%
%    lfp            list of (time,voltage1,...,voltageN) tuples
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

lfpstruct = getLFP(channels+1,varargin{:});
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
