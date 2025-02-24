function stats = MapStats2D_basic(map)

%MapStats - Compute statistics for a map Z = f(X,Y), or a curve Z = f(X).
%
%  Compute statistics on a continuous map, where one time-varying variable Z
%  is represented as a function of one or two time-varying variables X and Y.
%  The variable Z can either be a point process (typically, a list of spike
%  timestamps) or a continuous measure (e.g. the instantaneous velocity of
%  the animal, the spectral power of an LFP channel in a given frequency band,
%  the coherence between two oscillating LFP channels, etc.) Typical examples
%  of X and Y include spatial coordinates and angular directions.
%
%  USAGE
%
%    stats = MapStats2D_basic(map)
%
%    map            map obtained using <a href="matlab:help Map">Map</a>
%
%  OUTPUT
%
%    stats.specificity   spatial specificity (in bits, see Skaggs et al., 1993)
%    stats.informationPerSpike   spatial information, bits per spike; 
%    (Note that the specificity and informationPerSpike are very similar,
%    but have subtle difference. I am keeping both to make it consistent with
%    our previous calculations)
%    stats.informationPerSec   spatial information, bits per second
%    stats.sparsity   sparsity
%    stats.selectivity   peak rate/ mean rate ratio
%
% Original code from MapStats.m by Michaël Zugaro and PlceCellInfo.m by Azahara Oliva
% Modified by Wenbo Tang, Jan 31, 2023
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if ~isfield(map,'z'),
    map.z = map.rate;
    map = rmfield(map,'rate');
end

% Default values
stats.specificity = nan;
stats.informationPerSpike = nan;
stats.informationPerSec = nan;
stats.sparsity = nan;
stats.selectivity = nan;

if isempty(map.z), return; end

% Compute the spatial specificity of the map, based on the formula proposed by Skaggs et al. (1993):
%  specificity = SUM { p(i) . lambda(i)/lambda . log2(lambda(i)/lambda) }
% Compute the spatial information, sparsity and selectivity
T = sum(map.time(:));
p_i = map.time/(T+eps); % Probability of the animal occupying bin 'i'
lambda_i = map.z;
lambda = lambda_i(:)'*p_i(:);

meanFiringRate = sum(sum(map.z.*map.time))./T;
logArg = map.z./meanFiringRate;

if T == 0 || lambda == 0
    stats.specificity = nan; % set to nan instead of 0
    stats.informationPerSpike = nan;% set to nan instead of 0
    stats.informationPerSec = nan; 
    stats.sparsity = nan;
    stats.selectivity = nan;
else
    stats.specificity = sum(sum(p_i.*lambda_i/lambda.*log2(lambda_i/lambda)));
    logArg(logArg == 0) = 1;
    stats.informationPerSpike  = sum(sum(p_i.*logArg.*log2(logArg))); % bits per spike.
    stats.informationPerSec = sum(sum(p_i.*map.z.*log2(logArg))); % bits per second.
    stats.sparsity = ((sum(sum(p_i.*map.z))).^2)/sum(sum(p_i.*(map.z.^2)));
    stats.selectivity = max(max(map.z))./meanFiringRate;
end

