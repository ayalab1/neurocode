function [m,r,p] = meanangle(angles,weights,d)

% Computes the mean of a set of angles (in radians).
%
% Usage: [m,r,pvalue] = meanangle(angles,weights,dim/groups)
% where m is the mean angle, and r is the mean vector length
%
% 'angles' is a vector or matrix of angles (in radians).
% 'weights' is the relative weight given to each angle (default=ones)
% 'd' is either a dimension in the case of matrix input, or
% a grouping variable.
%
% EXAMPLE:
% [m,r,p] = meanangle(Phase(theta,spikes(:,1)),1,spikes(:,2));
% m,r, and p are vectors of n elements, where n is the largest spike id provided.
% EXAMPLE WITH WEIGHTS:
% [~,gAmplitude] = Phase(FilterLFP(lfp,'passband','gamma'));
% [m,r,p] = meanangle(Phase(FilterLFP(lfp,'passband','theta')),gAmplitude(:,2));
% will give the stats (m,r, and p) for gamma amplitude modulation to theta phase.

if nargin<3
    d = find(size(angles)>1,1);
    if isempty(d) %This is a scalar
        m = angles;
        return
    end
end


if ~exist('weights','var') || isempty(weights) || length(weights)==1,
    weights = ones(size(angles));
end

angles = exp(1i*angles).*weights;

if length(d)==1, % d = dimension
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

if nargout==2,
r = r.*(n-1)./n;
end




