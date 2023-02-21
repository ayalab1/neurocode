function bst = get_participation(st, events, par_type)
% get_participation computes the participation matrix for each cell given the 
% spike train and event timestamps.
%
% bst = get_participation(st, events) computes the participation matrix for 
% each cell given the spike train and event timestamps. The spike train is 
% a SpikeArray. The events argument is a IntervalArray.
%
% The function returns bst, which is an n_cells x n_events participation matrix.
% Each element of the matrix represents the number of spikes that occurred 
% for a particular cell during a particular event interval.
%
% bst = get_participation(st, events, par_type) computes the participation
% matrix for each cell given the spike train and event timestamps and the
% type of participation metric to calculate. The input argument `par_type` 
% specifies the type of participation metric to compute. It can be one of 
% three strings: 'counts' (default), 'binary', or 'firing_rate'. If `par_type`
% is set to 'counts', the participation matrix is simply the number of spikes
% for each cell during each event interval. If `par_type` is set to 'binary', 
% the participation matrix is a binary matrix where an element is 1 if there 
% is at least one spike during an event interval and 0 otherwise. If `par_type`
% is set to 'firing_rate', the participation matrix is the firing rate for 
% each cell during each event interval.
%
% Example:
%   % Generate SpikeArray and IntervalArray
%   st = SpikeArray(spikes.times)
%   ripples = IntervalArray(ripples.timestamps)
%   bst = get_participation(st, ripples, 'binary');
%   participation_prob = mean(bst,2)
% Inputs:
% - st: SpikeArray class.
% - events: IntervalArray class.
% - par_type: Type of participation metric to compute. Can be 'counts' 
%   (default), 'binary', or 'firing_rate'.
%
% Outputs:
% - bst: Participation matrix for each cell given the spike train and 
%   event timestamps.
%
% Ryan H 2023


if ~isequal(class(st),'SpikeArray')
   error('st needs to be SpikeArray class') 
end

if ~isequal(class(events),'IntervalArray')
   error('events needs to be IntervalArray class') 
end

if nargin < 3
    par_type = "counts";
end

% initialize the output matrix
bst = zeros(st.n_cells, events.n_intervals);

for i = 1:st.n_cells
    bst(i,:) = CountInIntervals(st.spikes(st.uid==i),events.intervals);
end

% calculate participation metric based on specified type
switch par_type
    case 'counts'
        % do nothing
    case 'binary'
        bst = (bst > 0);
    case 'firing_rate'
        bst = bst ./ events.lengths';
end
end