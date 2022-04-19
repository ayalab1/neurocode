function [indices,values] = FindClosest(reference, query, mode)

%FindClosest - field the indices of a reference vector closest to a query vector
% Look up a query vector in a reference vector, and for each value in 
% the query vector, return the index of the closest value in the
% reference vector. 
%
% USAGE
%
%    [indices,values] = FindClosest(reference, query, mode)
%
%    reference      reference signal which will be used as the table look-up
%    query          query signal. For each value of this signal, the function
%                   will look for the closest value in the reference signal
%    mode           a string ('either','higher', or 'lower') indicating 
%                   whether the function should find the closest value in
%                   reference higher than the query ('higher) or lower than
%                   the query ('lower') or the closest value regardless of
%                   the direction of the difference (default = 'either')
%
%
% EXAMPLE
% indices = FindClosests(spikes(:,1),deltas(:,2)) will return the indices of spikes closest to the delta wave peak
% d = deltas(:,2) - spikes(indices) will contain all the minimal distances between delta wave peaks and spikes (one distance per each delta wave)
%
% NOTE 
% This is an optimized implementation of computing the following trivial code:
% for i=1:length(reference), indices(i,1) = find(abs(query-reference(i)) == min(abs(query-reference(i))), 1, 'first'); end
% values = reference(indices);
%
% (C) 2019 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('mode','var'),mode = 'either'; end
if isempty(reference) || isempty(query), indices = []; values = []; return; end
if length(reference)==1,indices = 1; values=reference; return; end

nans = isnan(reference);
nonnans = find(~nans);
reference(nans) = [];

% Make sure values to be matched are unique
[u,i] = unique(reference);
% Find closest index
if strcmp(mode,'higher'),
    indices = ceil(interp1(u,(1:length(u))',query));
elseif strcmp(mode,'lower'),
    indices = floor(interp1(u,(1:length(u))',query));
else
    indices = round(interp1(u,(1:length(u))',query));
end

% Closest value to an undershooting is the smallest one in the reference list
indices(indices<1 | query<min(u))=1;
% Closest value to an overshooting query is the largest one in the reference list
indices(indices>length(i) | query>max(u)) = length(i);

indices = i(indices);
values = reference(indices);
indices = nonnans(indices);
