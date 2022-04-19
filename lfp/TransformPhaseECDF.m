function corrected = TransformPhaseECDF(phases,phases0)

%TransformPhaseECDF - correct phases for signal nonuniformity (e.g. asymmetry)
%
% TransformPhaseECDF corrects the 'phases' given to take into account
% the non-uniformity of the signal, as evidenced by 'phases0'
% Siapas 2005 correction for non-uniformity using the emprical cumulative distribution of phases
% Let FU :[−π,π)→[0,1) be the empirical cumulative distribution function corresponding 
% to ϕU. Then, Ψ(x) = 2πFU(x) − π yields the desired transformation. 
% This transformation is equivalent to the computation of the circular ranks or uniform 
% scores of the elements of ϕU. Clearly, if ϕU is uniform to begin with, Ψ reduces to the identity map.
%
% EXAMPLE USAGE:
% phases0 = Phase(filtered);
% phases = Phase(filtered,spikes(:,1));
% phases = TransformPhaseECDF(phases(:,2),phases0(:,2));
%
% Copyright (C) 2019 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure angles are between 0 and 2pi:
phases0 = wrap(phases0,2);
phases = wrap(phases,2);

% Get the empirical distribution
[ht,h] = ecdf(phases0);
[~,u,~] = unique(h);
corrected = interp1(h(u),ht(u),phases);
corrected = corrected * 2*pi;
