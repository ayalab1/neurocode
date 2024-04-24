function simulated = GeneratePoissonSpikes(spikes,maxTime)

% Simulate poisson spikes with the same firing rate at the provided spikes.
% Note: spikes should be provided in [timestamp id] format
%
%  INPUTS
%    spikes      firing rate of spikes to be simulated
%    maxTime     max time to simulate spikes over
%
%
% 
%  OUTPUTS
%    simulated    simulated poission spikes
%  SEE ALSO
%
%
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

simulated = [];
for i=1:max(spikes(:,2))
    s = spikes(spikes(:,2)==i);
    isi = exprnd(maxTime/length(s),length(s),1);
    sh = cumsum(isi);
    sh(:,2) = i;
    simulated = [simulated;sh];
end
simulated(simulated(:,1)>maxTime,:) = [];
simulated = sortrows(simulated);
