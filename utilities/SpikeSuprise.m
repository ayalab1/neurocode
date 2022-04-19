function [z,p] = SpikeSuprise(count, firingRates, binSize)
%SpikeSuprise - estimate the surprise of a spike count assuming poisson spiking
% Estimate probability of observing a number of spikes ('count')
% counted over a window of a given duration ('binSize')
% given a neuron's average firing rate ('firingRates')
% assuming Poisson spiking.
% The surprise is reported in z units where z>1.96
% corresponds to spike counts significantly higher than
% what would be expected from a Poisson neuron, and
% z<-1.96 corresponds to significantly lower spike
% counts than expected.
% 'count' and 'firingRates' can be vectors corresponding to
% multiple neurons
%
% Copyright (C) 2016 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

r = firingRates;
k = count;
dt = binSize;

% Code for p(x==k), later replaced by code for p(x>k)
% p = exp(k.*log(r.*dt)-logfactorial(k)-r.*dt); 
% smaller = k<r.*dt; % for the sign of z

lambda = r.*dt;
p = poisscdf(k,lambda);
smaller = p<0.5;

p(~smaller) = poisscdf(k(~smaller),lambda(~smaller),'upper');

z = p2z(p*2); % make it two-tailed. p=0.98 is just as unlikely as p=0.02;
z(abs(z)==Inf) = 40; % matlab rounds p-values smaller than z=40 as 0, i.e. z=Inf
z(smaller) = -z(smaller);
