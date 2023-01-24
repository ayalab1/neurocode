function [spectra, f, ns] = RelativePointSpectra(spikes,cycles)

% Compute point spectra relative to an oscillation
% This is useful for computing phase precession of a cell relative to the ongoing theta oscillation
%
% EXAMPLE
% [spectra, f, ns] = RelativePointSpectra(cell_metrics.spikes.times,Restrict(thetacycles.timestamps,behavior.trials));
% PlotColorMap(nanzscore(Smooth(spectra,[0 2]),[],2),'x',f); PlotHVLines(1,'v','k'); xlabel('frequency (relative to cycles)'); ylabel('unit ID'); set(get(colorbar,'YLabel'),'String','power (z units');
% Copyright (C) 2022-2023 by Ralitsa Todorova
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% make sure spikes are a matrix:
if iscell(spikes), spikes =  Group(spikes{:},'sort',true); end
id = spikes(:,2);
nUnits = max(id);

f = linspace(0.5,1.5,1000);
cycletime = RelativeTimes(spikes(:,1),cycles,[0 1]); % get "theta phase"
[~,w] = InIntervals(spikes(:,1),cycles);
cycletime = cycletime+w-1;
ok = ~isnan(cycletime);
timesCell = cell(nUnits,1); for i=1:nUnits, timesCell{i} = cycletime(id==i & ok); end
[spectra0,f0] = PointSpectra(timesCell,'range',[0 3],'pad',0,'frequency',100);
spectra = interp1(f0,spectra0,f)';
ns = cellfun(@length,timesCell)'; % number of spikes used to compute the spectrogram (useful if the user wants to remove units of few spikes post-hoc)
