function [d,s,linkage] = FunctionalClustering(spikes,varargin)

%FunctionalClustering - Determine spike train similarity tree using FCA analysis.
%
% Determine which spike trains are most similar in a neuronal population. This
% uses the Functional Clustering Algorithm (FCA) of Feldt et al. (2009).
%
% For each pair of neurons i and j, we first measure the distance (time elapsed)
% between each spike of i and the closest spike of j, and use the average of these
% distances as an estimate of the similarity between i and j. We then swap the
% two neurons and measure the similarity between j and i. The average minimal
% distance (AMD) is defined as the mean between the two similarity values. Here,
% we use a normalizing factor so that a pair of independent Poisson processes
% has an AMD of 1.
%
% To define the scaled significance, the spike trains are repeatedly jittered by
% a random amount, and the distribution of jittered spike trains is determined.
% The AMD is centered and normalized using the median and 95% percentile of this
% distribution. Thus, a scaled significance of 1 corresponds to a p value of 0.05.
%
% Finally, the similarity tree is computed recursively by finding the two neurons
% with the highest AMD, grouping them into a new cluster, and repeating until all
% neurons have been grouped together.
%
%  USAGE
%
%    [d,s,linkage] = FunctionalClustering(spikes,<options>)
%
%    spikes         two-column matrix of timestamps and unit IDs,
%                   provided by <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a> using 'output' = 'numbered'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'jitter'      jitter standard deviation, in s (default = 5)
%     'nJitters'    number of jittered spike trains (default = 1000)
%    =========================================================================
%
%  OUTPUT
%
%    d              normalized average minimal distance for each unit pair
%    s              scaled significance for each unit pair
%    linkage        network connectivity tree
%
%  SEE
%
%  See also PlotLinkage.

% Copyright (C) 2015 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
jitterSD = 5;
nJitters = 1000;

% Check parameter
if ~isdmatrix(spikes,'@2') || ~isivector(spikes(:,2)),
	error('Incorrect spikes (type ''help <a href="matlab:help FunctionalClustering">FunctionalClustering</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FunctionalClustering">FunctionalClustering</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'jitter',
			jitterSD = varargin{i+1};
			if ~isdscalar(jitterSD,'>0'),
				error('Incorrect value for property ''jitter'' (type ''help <a href="matlab:help FunctionalClustering">FunctionalClustering</a>'' for details).');
			end
		case 'njitters',
			nJitters = varargin{i+1};
			if ~isiscalar(nJitters,'>0'),
				error('Incorrect value for property ''nJitters'' (type ''help <a href="matlab:help FunctionalClustering">FunctionalClustering</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionalClustering">FunctionalClustering</a>'' for details).']);
	end
end

spikes = sortrows(spikes);
t = spikes(:,1);
id = spikes(:,2);
n = max(id);
T = max(t) - min(t);
% Make sure jittered spikes remain within the time window if this is small
jitterSD = min([jitterSD T/10]);

% The following two matrices have sizes 2n x 2n (rather than n x n) because this will speed up
% computations in the second part of the algorithm, but only the upper n x n square is actually
% used to compute the AMD and scaled significance of the unit pairs
AMD = nan(2*n);
scaledSignificance = nan(2*n);
linkage = nan(1);
jitteredAMD = nan(n,n,nJitters);

% Compute average firing rate for each unit
for i = 1:n,
	f(i,1) = sum(id==i)/T;
end

% Compute normalized average minimal distances for each pair of units
for i = 1:n,
	for j = i+1:n,
		AMD(i,j) = ComputeAMD(spikes,f,i,j);
	end
end
d = AMD(1:n,1:n);

% Compute normalized average minimal distances for each pair of jittered units (for each jitter)
for jitter = 1:nJitters,
	jittered = sortrows([t+jitterSD*randn(size(t)) id]);
	for i = 1:n,
		for j = i+1:n,
			jitteredAMD(i,j,jitter) = ComputeAMD(jittered,f,i,j);
		end
	end
end

% Center and reduce to get the scaled significance
medians = quantile(jitteredAMD,0.5,3);
scale = quantile(jitteredAMD,0.05,3) - medians;
s = (AMD(1:n,1:n)-medians)./scale;

% Compute linkage only if required (time consuming)
if nargout < 3, return; end

% Merge most similar spike trains and update matrix
scaledSignificance(1:n,1:n) = s;
linkage = nan(numel(unique(id))-1,4);
iteration = 0;
jitteredAMD = nan(n,nJitters);
while numel(unique(id))>1,
	iteration = iteration+1;
	% Find max scaled significance
	ss = max(scaledSignificance(:));
	[i,j] = find(scaledSignificance==ss,1,'first');
	% Drop the corresponding line and column
	scaledSignificance(:,[i j]) = nan;
	scaledSignificance([i j],:) = nan;
	% Update connectivity tree
	linkage(iteration,:) = [i j ss AMD(i,j)];
	% Merge these units into a new cluster
	n = n+1;
	spikes(id==i|id==j,2) = n;
	id = spikes(:,2);
	f(n,1) = sum(id==n)/T;
	
	% update AMD and scaled significance
	for jitter = 1:nJitters,
		jittered = sortrows([t+jitterSD*randn(size(t)) id]);
		for i = 1:n-1,
			jitteredAMD(i,jitter) = ComputeAMD(jittered,f,i,n);
		end
	end
	for i=1:n-1,
		AMD(i,n) = ComputeAMD(spikes,f,i,n);;
		medians = quantile(jitteredAMD(i,:),0.5,2);
		scale = quantile(jitteredAMD(i,:),0.05,2) - medians;
		scaledSignificance(i,n) = (AMD(i,n)-medians)./scale;
	end
end

% ------------------------------- Helper function -------------------------------

function AMD = ComputeAMD(spikes,f,i,j)

AMD = nan;
spikes_ij = spikes(ismember(spikes(:,2),[i j]),:);
if numel(unique(spikes_ij(:,2))) < 2, return; end;
% Distances between successive spikes, whatever unit they belong to
D = [1000 1000; diff(spikes_ij); 1000 1000];
distance = D(:,1);
different = logical(D(:,2)); % same or different unit?
% Distance to previous spike of the other unit
previous = CumSum(distance,different);
% Distance to next spike of the other unit (flip vectors upside down)
next = flipud(CumSum(flipud(distance),flipud(different)));
% Minimum between the two (for each spike)
m = min([[NaN; previous(2:end-1)] [next(2:end-1); NaN]],[],2);
% (Normalized) average minimal distance (mean of average minimal distances from i to j and from j to i)
d_ij = mean(m(spikes_ij(:,2)==i));
d_ji = mean(m(spikes_ij(:,2)==j));
AMD = mean([d_ij d_ji]) / ((f(i)+f(j))/(4*f(i)*f(j)));
