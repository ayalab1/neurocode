function [consolidated,target] = ConsolidateIntervalsFast(intervals,varargin)

%ConsolidateIntervalsFast - a faster version of ConsolidateIntervals.
% Attention - this function assumes the intervals are sorted in time and that
% while intervals may overlap, no interval is completely contained within another 
% interval. If these assumptions are not met, use ConsolidateIntervals instead.
%
% Consolidate overlapping intervals, e.g. replace [10,20] [15,25] with [10,25].
%
%  USAGE
%
%    [consolidated,target] = ConsolidateIntervals(intervals,<options>)
%
%    intervals      list of intervals
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'strict'      intervals with common bounds are consolidated ('off')
%                   or kept separate ('on') (default = 'off')
%     'epsilon'     intervals with close enough bounds (distance lesser than
%                   epsilon) are also consolidated (default = 0)
%    =========================================================================
%
%  OUTPUT
%
%    consolidated   consolidated intervals
%    target         for each original interval, the index of the consolidated
%                   interval to which it belongs (empty intervals yield NaN)
%
%  SEE
%
%    See also SubtractIntervals, ExcludeIntervals, InIntervals, Restrict,
%    FindInInterval, CountInIntervals, PlotIntervals.


% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
epsilon = 0;
strict = 'off';

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'strict',
			strict = lower(varargin{i+1});
			if ~isstring(strict,'on','off'),
				error('Incorrect value for property ''strict'' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
			end
		case 'epsilon',
			epsilon = varargin{i+1};
			if ~isdscalar(epsilon,'>0'),
				error('Incorrect value for property ''epsilon'' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).']);
	end
end

% Get difference between two consequitive intervals
vector = reshape(intervals',[],1);
d = diff(vector);
d = d(2:2:end);

% Define problematic spots
if strcmp(strict,'on'),problematic = d<epsilon;else problematic = d<=epsilon; end

% Leave function if intervals don't need consolidating
if ~any(problematic),
    consolidated = intervals;
    target = (1:size(intervals,1))';
    return
end

% Mark already consolidated intervals
done = intervals([(~problematic);true] & [true;(~problematic)],:);

% For problematic regions, use start of the first interval and stop of the last
indices = ToIntervals(problematic); 
indices(:,2) = indices(:,2)+1;
fixed = [intervals(indices(:,1),1) intervals(indices(:,2),2)];

% Put them back together
consolidated = sortrows([done; fixed]);

% If asked, provide indices of new intervals in which old intervals fall
if nargout<1,
    [~,target] = InIntervals(mean(intervals,2),consolidated); % using interval midpoint as a proxy for the whole interval
end