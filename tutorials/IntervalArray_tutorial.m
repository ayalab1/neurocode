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

disp("total duration: ");
ripple_epochs.duration

disp("minimum bound of all intervals: ");
ripple_epochs.min()

disp("maximum bound of all intervals: ");
ripple_epochs.max()

% you can return a copy of the class at a particular index
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
ripple_epochs = ripple_epochs.intersect(beh_epochs(1))

%% check that these ripples are within the interval

ripple_epochs.max <= beh_epochs(1).max

ripple_epochs.min >= beh_epochs(1).min

%% restrict to nrem within epoch 1

nrem_epochs = IntervalArray(SleepState.ints.NREMstate);

ripple_epochs = ripple_epochs.intersect(nrem_epochs)


