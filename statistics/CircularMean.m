function [m,r,p] = CircularMean(angles,weights,d)

% CircularMean - Computes the mean of sets of angles (in radians).
%
%
%  USAGE
%
%    [m,r,pvalue] = CircularMean(angles,weights,dim/groups)
%
%    angles         vector or matrix of angles (in radians)
%    weights        relative weight given to each angle (default = ones)
%    dim            optional dimension along which the mean should be computed 
%                   (required if input is matrices rather than vectors)
%    groups         optional grouping variable. The results for each group
%                   will be computed separately
%
%  OUTPUT
%
%    m              mean angle
%    r              resultant vector length
%    p              p-value (Rayleigh test)
%
%
%  EXAMPLE
%
% [m,r,p] = CircularMean(Phase(theta,spikes(:,1)),1,spikes(:,2));
% m,r, and p are vectors of n elements, where n is the largest spike id provided.
%
% [~,gAmplitude] = Phase(FilterLFP(lfp,'passband','gamma'));
% phases = Phase(FilterLFP(lfp,'passband','theta'))
% [m,r,p] = CircularMean(phases(:,2),gAmplitude(:,2));
% will give the stats (m,r, and p) for gamma amplitude modulation to theta phase.
%
% Copyright (C) 2004-2022 by Michaël Zugaro & Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin<3
    d = find(size(angles)>1,1);
    if isempty(d) %This is a scalar
        m = angles;
        return
    end
end


if ~exist('weights','var') || isempty(weights) || length(weights)==1
    weights = ones(size(angles));
end

angles = exp(1i*angles).*weights;

if length(d)==1 % d = dimension
    mid = nanmean(angles,d);
    n = size(angles,d);
else
    % d = group
    n = accumarray(d,1);
    mid = accumarray(d,angles)./n;
end

m = atan2(imag(mid),real(mid));
r = sqrt(imag(mid).^2+real(mid).^2);
p = exp(sqrt(1+4.*n+4*(n.^2-(n.*r).^2))-(1+2.*n)); % Zar, Biostatistical Analysis, p. 617
end



