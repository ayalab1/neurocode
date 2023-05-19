function H = RasterPlot(spikes, height, varargin)

%RasterPlot - create raster plot from spike times and IDs
%
% "spikes" is a list of [timestamps id]
% provide desired height of spikes
% any other inputs will be passed on to "plot"
%
% Copyright (C) 2018-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin==1,
    height = 1;
end
    
times = spikes(:,1);
rows = spikes(:,2);

times = [times times nan(size(times))]';
rows =  [rows-0.45*height rows+0.45*height nan(size(rows))]';

H = figure();
plot(times(:),rows(:),varargin{:});