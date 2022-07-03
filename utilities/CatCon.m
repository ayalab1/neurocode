function [start, stop, dur, peak, numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenation with Condition function
% This function is to be used in conjunction with the find_HSE_b code. For
% proper use, it will intake event start, stop, and peak times, as well as
% the amplitude of each peak and a logical array designating what must be
% concatenated. This must be calculated outside of the function. The fields 
% that aren't needed can be bypassed by inputting arrays of zeros/ones.
%
%%% INPUTS %%%
% evtstart:     Array with the start times of the events of interest
% evtstop:      Array with the stop times of the events of interest
% evtpeak:      Array with the peak times of the events of interest
% evtamp:       Array with the amplitudes of the peak for each event
% flagConc:     Logical array denoting whether values should be
%               concatenated. Note that it must be arranged so that [1 0]
%               designates that event 1 and 2 should be concatenated, but 2
%               and 3 should not. flagConc will have a length of
%               length(evtstart)-1.
%
%%% OUTPUTS %%%
% start:        An updated array of start times of the events after
%               concatenation
% stop:         An updated array of stop times of the events after
%               concatenation
% dur:          An updated array of event durations after concatenation
% peak:         An updated array of peak times of the events after
%               concatenation
% numCat:       Array containing the number of events concatenated within
%               each new event. Useful for tracking burst length if
%               concatenating ISIs. 
%
% Lindsay Karaba, 2021, AYA Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

concHSE_start = [];
concHSE_stop = [];
concHSE_peak = [];
concHSE_count = [];
cHSE = 1;
catting = 0; %are we in the process of concatenating - to start, no
catCnt = 1;
    
for i = 1:length(evtstart)-1 % iterate through start times
    if ~catting %if we are not already concatenating, we need to set our next start time
        concHSE_start(cHSE) = evtstart(i);
        tempPeak = [];
        tempPeak(1,1) = evtpeak(i);
        tempPeak(2,1) = evtamp(i);
        catCnt = 1;
    end
    %check if we should be concatenating or if we should stop
    if flagConc(i) %are we flagging the next set to keep 
        catting = 1;
        tempPeak(1, end+1) = evtpeak(i); %we might double count but that's okay
        tempPeak(2, end) = evtamp(i);
        catCnt = catCnt + 1;
    else
        catting = 0;
        tempPeak(1, end+1) = evtpeak(i); %we might double count but that's okay
        tempPeak(2, end) = evtamp(i);
        concHSE_stop(cHSE) = evtstop(i);
        maxPeak = [];
        maxPeak = find(tempPeak(2,:) == max(tempPeak(2,:)));
        concHSE_peak(cHSE) = tempPeak(1,maxPeak(1));
        concHSE_count(cHSE) = catCnt;
        cHSE = cHSE+1;
    end
end

if length(concHSE_stop) < length(concHSE_start) 
    tempPeak(1, end+1) = evtpeak(i); %we might double count but that's okay
    tempPeak(2, end) = evtamp(i);
    concHSE_stop(cHSE) = evtstop(i);
    maxPeak = [];
    maxPeak = find(tempPeak(2,:) == max(tempPeak(2,:)));
    concHSE_peak(cHSE) = tempPeak(1,maxPeak(1));
    concHSE_count(cHSE) = catCnt;
    cHSE = cHSE+1;
end

start = concHSE_start;
stop = concHSE_stop;
dur = stop-start;
peak = concHSE_peak;
numCat = concHSE_count;

end