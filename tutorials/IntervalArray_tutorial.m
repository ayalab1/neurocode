% Demo of IntervalArray
%
% IntervalArray class represents an array of intervals
% 
% Provides various methods for manipulating and analyzing intervals, 
%   such as expanding, shrinking, finding intersections, and unions
% 
% Can validate intervals for errors and sort them by start time
% 
% Can return various properties of intervals, such as start and stop times,
%   lengths, and total duration
% 
% Can remove empty intervals, check for equality, and check if a point is 
%   within one of the intervals
% 
% Can return a copy of the class with just one interval using the () operator.
% 
% You can also use MATLAB Operators! 
%   & for intersection
%   | for union
%   - for set difference
%   + for set addition
%   ~ for complement
%   
% Can be useful in a variety of applications such as signal processing 
%   where intervals are frequently used.

%%
% load in ripples, session, sleepstates
load('Z:\Data\HMC2\day9\day9.ripples.events.mat')
load('Z:\Data\HMC2\day9\day9.session.mat')
load('Z:\Data\HMC2\day9\day9.SleepState.states.mat')

%% add ripples to IntervalArray class
ripple_epochs = IntervalArray(ripples.timestamps)

%% demo types of simple methods
disp("number of intervals: ");
ripple_epochs.n_intervals

% check the overall duration (all intervals added up)
disp("total duration: ");
ripple_epochs.duration

% also printing out the class will produce a nice overview of data
disp('look at pretty disp')
ripple_epochs

% first start
fprintf('\n')
disp("minimum bound of all intervals: ");
ripple_epochs.min()

% last stop
disp("maximum bound of all intervals: ");
ripple_epochs.max()

% you can easily return a copy of the class at a particular index
disp("retuning class by index: ");
ripple_epochs(2)

%%
beh_epochs = [];
for ep_i = 1:length(session.epochs)
    beh_epochs(ep_i,:) = [session.epochs{ep_i}.startTime,...
        session.epochs{ep_i}.stopTime];
end
beh_epochs = IntervalArray(beh_epochs);

% display plot using class method
beh_epochs.plot();

%% get the ripples within the first epoch
ripple_epochs_new = ripple_epochs.intersect(beh_epochs(1))

%% check that these ripples are within the interval

ripple_epochs_new.max <= beh_epochs(1).max

ripple_epochs_new.min >= beh_epochs(1).min

%% restrict to nrem

nrem_epochs = IntervalArray(SleepState.ints.NREMstate);

ripple_epochs_new = ripple_epochs.intersect(nrem_epochs)

%% you can use MATLAB Operators (& | + - ~) as well! 

% get intersection of multiple intervals using &
intervals_1 = ripple_epochs & nrem_epochs & beh_epochs(1)

% gives the same as this
intervals_2 = ripple_epochs.intersect(nrem_epochs).intersect(beh_epochs(1))

% lets check with ==
intervals_1 == intervals_2

%% Easily get complement with .complement or ~

disp('all epochs outside ripples')
~ripple_epochs

%% get union with .union or |
ripple_epochs | nrem_epochs

%% addition wth .plus or +
ripple_epochs + nrem_epochs

%% set difference with .setdiff or -
ripple_epochs - nrem_epochs


