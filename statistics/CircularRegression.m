function [beta,R2,p] = CircularRegression(x,angles,p,varargin)

%CircularRegression - Linear-circular regression.
%
%  There are three possible methods to estimate the regression parameters:
%
%  1) Use the method of Kempter et al. (2012). The slope maximizes the resultant length
%     of the residual distribution, the phase minimizes the circular variance of the
%     residual distribution.
%
%  2) Fit 'barber-pole' lines to the data. The model consists of three parallel lines
%     y = a.x + b + k*2pi (k in {-1,0,1}). The distance of a point to the model is
%     defined as the minimal distance to either of the three lines. The best-fit model
%     is computed using a least squared error approach.
%
%  3) Compute Thiel-Sen estimators after reorganizing angles around the regression line
%     obtained above (move by multiples of 2.pi).
%
%  Regression can be performed either on a list of (x,phi) pairs, or using a probability
%  distribution of phi for each value of x.
%
%  USAGE
%
%    [beta,R2,p] = CircularRegression(x,angles,<options>)
%
%    x              values for linear independent variable
%    angles         values for circular dependent variable (in radians)
%    <options>      optional list of property-value pairs (see table below)
%
%    [beta,R2,p] = CircularRegression(x,angles,p,<options>)
%
%    x              N linear variable bins
%    angles         M circular variable bins
%    p              MxN matrix where the nth column is the angular probability
%                   distribution for the nth value of x (angular bins appear
%                   in ascending order from top to bottom)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'slope'       search start value for slope (default = 0)
%     'method'      regression method: 'k' for Kempter (default), 'pole' for
%                   barber pole, or 'ts' for Theil-Sen (but 'ts' cannot be
%                   computed for probability distribution matrices)
%    =========================================================================
%
%  OUTPUT
%
%    beta           regression slope and intercept
%    R2             coefficient of determination
%    p              p-value associated with H0: slope = 0 (Kempter only)
%
%  SEE
%
%    See also CircularVariance, CircularConfidenceIntervals, Concentration,
%    ConcentrationTest.

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
method = 'k';
slope = 0;

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end

if ~isdvector(x),
	error('Incorrect x values (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end
if ~isdvector(angles),
	error('Incorrect angles (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end

if nargin >= 3
	if isdmatrix(p),
		% Two lists of bins and a probability matrix
		if size(p,1) ~= length(angles) || size(p,2) ~= length(x),
			error('Incoherent parameter sizes (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
		end
	else
		% Two lists of samples (x and angles)
		varargin = {p,varargin{:}};
		p = [];
	end
else
	p = [];
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'method',
			method = varargin{i+1};
			if ~isastring(method,'pole','k','ts'),
				error('Incorrect value for property ''method'' (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
			end
		case 'slope',
			slope = varargin{i+1};
			if ~isdscalar(slope),
				error('Incorrect value for property ''slope'' (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
			end
		otherwise,
			parameters = varargin{i:end};
			if ~isa(parameters,'cell'), parameters = {parameters}; end
			break;
	end
end

if strcmp(method,'ts') && ~isempty(p),
	error('Thiel-Sen regression cannot be computed for probability distribution matrices (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end

isradians(angles);
if ~isempty(p),
	% Replicate x values (one copy for each angular bin)
	[m,n] = size(p);
	x = reshape(x,1,length(x));
	x = repmat(x,m,1);
	x = x(:);
	% Angular bins (one copy for each x bin)
	phi = angles(:);
	phi = repmat(phi,n,1);
	prob = p(:);
else
	x = x(:);
	phi = angles(:);
	prob = ones(size(phi));
	prob = prob/length(prob);
end

beta = [nan nan];
R2 = nan;
p = nan;

if isempty(angles) || isempty(x), return; end

% Store data in global variables (necessary for minimization)
global CircularRegression_x CircularRegression_phi CircularRegression_prob;
CircularRegression_x = x;
CircularRegression_phi = phi;
CircularRegression_prob = prob;

% ------------------------------------------------------------------------------------------------
% Compute barber-pole regression
% (this is a first step for Thiel-Sen)

if strcmp(method,'pole') || strcmp(method,'ts'),
	% Minimize residual error (least squared error)
	[beta,RSE] = fminsearch(@ResidualSquaredError,[slope 0]);

	% Compute coefficient of determination
	TSE = norm(prob.*(phi-CircularMean(phi)))^2;
	R2 = 1-RSE/TSE;
	% p-value not implemented
end

% ------------------------------------------------------------------------------------------------
% Compute 'Kempter' regression

if strcmp(method,'k'),
	% Maximize resultant length
	[a,R] = fminsearch(@ResultantLength,slope);

	% Estimate phase shift
	C = sum(prob.*cos(phi-a*x));
	S = sum(prob.*sin(phi-a*x));
	phi0 = atan2(S,C);

	beta = [a phi0];

	% Compute coefficient of determination
	theta = mod(abs(a)*x,2*pi);
	phi_bar = atan2(sum(prob.*sin(phi)),sum(prob.*cos(phi)));
	theta_bar = atan2(sum(prob.*sin(theta)),sum(prob.*cos(theta)));
	N = sum(prob.*sin(phi-phi_bar).*sin(theta-theta_bar));
	D = sqrt(sum(prob.*sin(phi-phi_bar).^2)*sum(prob.*sin(theta-theta_bar).^2));
	rho = N/D;
	R2 = rho^2;

	% Compute p-value
	n = length(phi);
	lambda_02 = 1/n*sum(prob.*sin(theta-theta_bar).^2);
	lambda_20 = 1/n*sum(prob.*sin(phi-phi_bar).^2);
	lambda_22 = 1/n*sum(prob.*sin(phi-phi_bar).^2.*prob.*sin(theta-theta_bar).^2);
	z = rho*sqrt(n*lambda_20*lambda_02/lambda_22);
	p = 1-erf(abs(z)/sqrt(2));
	
end

% ------------------------------------------------------------------------------------------------
% Compute Theil-Sen estimator

if strcmp(method,'ts'),
	% Move phi values by multiples of 2.pi to get them as close as possible to the regression line
	d = phi-beta(1)*x+beta(2);
	rd = sign(d).*floor(abs(d/(2*pi)));
	phi = phi+rd*2*pi;

	% Keep at most 500 random samples
	n = length(x);
	if n == 1,
		beta = [nan nan];
		R2 = nan;
		return;
	end
	pp = randperm(n);
	n = min([n 500]);
	x0 = x(pp(1:n));
	phi0 = phi(pp(1:n));

	% Compute pairwise slopes
	slopes = nan(nchoosek(n,2),1);
	k = 1;
	for i = 1:n-1,
		for j = (i+1):n,
			slopes(k) = (phi0(i)-phi0(j))/(x0(i)-x0(j));
			k = k + 1;
		end
	end

	% Take median
	beta(1) = median(slopes);
	beta(2) = median(phi-beta(1)*x);
	% Compute coefficient of determination
	RSE = sum((phi-(beta(1)*x+beta(2))).^2);
	R2 = 1-RSE/TSE;
	% p-value not implemented
end

% --------------------------------------------------------------------------------------------------------------

function RSE = ResidualSquaredError(beta)

% Retrieve data
global CircularRegression_x CircularRegression_phi;
x = CircularRegression_x;
phi = CircularRegression_phi;

% Model parameters: slope and intercept
a = beta(1);
b = beta(2);

% Three lines
y1 = a*x+b;
y2 = y1+2*pi;
y3 = y1-2*pi;

% Squared error
d = min(prob.*[(phi-y1).^2 (phi-y2).^2 (phi-y3).^2],[],2);
RSE = sum(d);

% --------------------------------------------------------------------------------------------------------------

function R = ResultantLength(beta)

% Retrieve data
global CircularRegression_x CircularRegression_phi CircularRegression_prob;
x = CircularRegression_x;
phi = CircularRegression_phi;
prob = CircularRegression_prob;

% Model parameters: slope
a = beta(1);

n = length(x);
C = 1/n*sum(prob.*cos(phi-a*x));
S = 1/n*sum(prob.*sin(phi-a*x));
R = sqrt(C.^2+S.^2);

% We must actually return -R because Matlab can only minimize (but not maximize...)
R = -R;
