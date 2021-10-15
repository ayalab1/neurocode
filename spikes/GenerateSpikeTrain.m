function [spikes,V] = GenerateSpikeTrain(I,threshold,varargin);

%GenerateSpikeTrain - Generate spike train using a leaky integrate and fire neuron model.
%
%  USAGE
%
%    [spikes,V] = GenerateSpikeTrain(I,threshold,<options>)
%
%    I              input current <a href="matlab:help samples">samples</a>
%    threshold      spiking threshold
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'delay'       delay in seconds (default = 0)
%     'tau'         decay rate (default = 0.015)
%     'Vr'          resting membrane potential (default = 0)
%     'sigma'       noise standard deviation (default = 0, no noise)
%     'R'           input resistance
%    =========================================================================
%
%  OUTPUT
%
%    spikes         list of spike timestamps
%    V              membrane potentials across time
%
%
%  EXAMPLE
%
%    To get spike trains of two neurons receiving correlated inputs, provide
%      I = mu + sigma(sqrt(1-c)*E + sqrt(c)*C)
%    where c is the correlation, E is the exclusive input (or noise) and C is
%    the common signal.


% Copyright (C) 2015-2016 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
delay = 0;
Vr = 0;
tau = 0.015;
sigma = 0;
R = 1;

% Check parameters
if ~isdmatrix(I,'@2'),
	error('Incorrect input current (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
end
if ~isdscalar(threshold),
	error('Incorrect spike threshold (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
end

% Browse parameters
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Optional parameter ' num2str(i) ' is not a valid property (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'delay',
			delay = varargin{i+1};
			if ~isdscalar(delay,'>=0'),
				error('Incorrect value for property ''delay'' (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
			end
		case 'tau',
			tau = varargin{i+1};
			if ~isdscalar(tau,'>=0'),
				error('Incorrect value for property ''tau'' (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
			end
		case 'vr',
			Vr = varargin{i+1};
			if ~isdscalar(Vr),
				error('Incorrect value for property ''Vr'' (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
			end
		case 'sigma',
			sigma = varargin{i+1};
			if ~isdscalar(sigma,'>=0'),
				error('Incorrect value for property ''sigma'' (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
			end
		case 'r',
			R = varargin{i+1};
			if ~isdscalar(R),
				error('Incorrect value for property ''R'' (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GenerateSpikeTrain">GenerateSpikeTrain</a>'' for details).']);
	end
end

% Initialization
Vmin = -threshold;
if tau == 0,
	warning('The value of tau cannot be exactly 0 (division by zero) - using ''realmin'' instead');
	tau = realmin;
end
dt = median(diff(I(:,1)));
m = min(I(:,1));
M = max(I(:,1));
t = (m:dt:M+dt/100)'; % +dt/100 added for underflow errors
delay = round(delay/dt);

% Create Gaussian noise
randn('seed',sum(100*clock)); % set seed for random number generator
n = length(t);
noise = sigma*sqrt(1/dt)*randn(n,1);

% Make sure I is evenly sampled (interpolate)
I = Interpolate(I,t);

% Start
spikes = [];
V = nan(n,1);
V(1) = Vr;

for i = 1:n-1,
    if V(i) > threshold,
        spikes = [spikes;t(i)];
        V(i) = Vr;
    elseif V(i) < Vmin,
        V(i) = Vmin;
    end
    if i+delay > length(I),
        V(i+1) = V(i) + dt*(-(V(i) - Vr)/tau);
    else
        V(i+1) = V(i) + dt*((R*I(i+delay,2)-(V(i) - Vr))/tau + noise(i+delay)/sqrt(tau));
    end
end
