function [m,boundaries] = CircularConfidenceIntervals(angles,alpha,nBootstrap)

%CircularConfidenceIntervals - Compute circular mean and confidence intervals.
%
% Compute circular mean and confidence intervals of circular data (see Fisher,
% Analysis of circular data, p. 88). If fewer than 25 angle measurements are
% provided, bootstrap is used.
%
%  USAGE
%
%    [mean,boundaries] = CircularConfidenceIntervals(angles,alpha,nBootstrap)
% Input
%    angles         angles in radians
%    alpha          optional confidence level (default = 0.05)
%    nBootstrap     optional number of bootstraps (for < 25 angle values)
%                   (default = 2000)
% Output
%   m               circular mean
%   boundaries      circular confidence intervals
%
%  SEE
%
%    See also CircularMean, CircularVariance, Concentration, ConcentrationTest.
%
%  Dependencies: CircularMean

% Copyright (C) 2004-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1
	error('Incorrect number of parameters (type ''help <a href="matlab:help CircularConfidenceIntervals">CircularConfidenceIntervals</a>'' for details).');
end
isradians(angles);
if isdscalar(angles)
	m = angles;
	boundaries = [m m];
	return
end

m = [];
boundaries = [];
n = length(angles);
if n == 0, return; end

if nargin < 2
	alpha = 0.05;
end
if nargin < 3
	if n < 25
		nBootstrap = 2000;
	else
		nBootstrap = 0;
	end
end

if nBootstrap == 0
	[m,r1] = CircularMean(angles);
	[~,r2] = CircularMean(wrap(2*angles));
	delta = (1-r2)./(2*r1.^2);
	sigma = sqrt(delta/n);
	sinarg = norminv(1-alpha/2)*sigma;
	if sinarg < 1
		err = asin(sinarg);
	else
		err = pi;
	end
	boundaries = [m-err; m+err];
else
	m = CircularMean(angles);
	b = bootstrp(nBootstrap,'CircularMean',angles);

	% unwrap data around mean value
	unwrapped = mod(b-repmat(m,nBootstrap,1)+pi,2*pi)+repmat(m,nBootstrap,1)-pi;

	boundaries(1,:) = prctile(unwrapped,100*(alpha/2));
	boundaries(2,:) = prctile(unwrapped,100-100*(alpha/2));
end