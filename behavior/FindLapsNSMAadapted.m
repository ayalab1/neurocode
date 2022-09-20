function laps=FindLapsNSMAadapted(Vts,Vdata,newLapThreshold)
%
% [laps, Vhs] = FindLaps_HorseShoe(V,newLapThreshold);
%
% Find Laps in "Horseshoe-Geometry" of a circular track
%
% INPUT:
% V    ...  tsd of angular position in [0,360] degree range, with horsehoe
%           opening at 0 degrees.
%
% newLapThreshold ... OPTIONAL endpoint proximity threshold in percent of track length (default = 15%);
%                     whenever rat enters the proximity zone of e.g. 15% of tracklength near a horseshoe end, a new lap
%                     is started and the maximum (or minimum) is searched
%                     for a Lap-Top (around 360 end)  or Lap-Bottom (around 0 end).
%
% OUTPUT:
% laps  .... 1*nLaps struct array with fields
%   laps(i).start_ts  ... start timestamp of i-th lap (in 0.1 millisec
%                              NSMA units)
%   laps(i).pos       ... the value of input position V (in circle geometry not horseshoe!) at lap start point (in degrees)
%   laps(i).start_idx ... the index of the new lap start frame in input V tsd
%   laps(i).direction ... +1/-1 for up/down laps (with increasing/decreasing position in angle degrees)
%
%  Vhs ....  tsd nFrames*1 vector of position in horseshoe geometry, where
%               negative direction laps are ranging from [-360,0] and positive direction
%               laps are ranging from [0 360]. NOTE: there are discontinuous jumps at
%               the horseshoe turn points around 0 degrees and around
%               360/-360 degrees. NOTE2: nFrames = length(Data(V))
%
%
%   Author: PL
%   Version: 0.9  05/12/2005
%   edited by Ryan Harvey to work with standard linear track

if nargin == 1
    newLapThreshold = 15;   % default value
end

% Vdata = data_video(:,2);
% Vts = data_video(:,1);

TL = abs(max(Vdata) - min(Vdata));   % track length in degrees
th1= min(Vdata)+TL*newLapThreshold/100;        % lower threshold for lower horeshoe end (degrees)
th2 = max(Vdata)-TL*newLapThreshold/100;      % upper threshold for upper horseshoe end (degrees)


% loop over all frames
laps(1).start_ts = Vts(1);
laps(1).pos = Vdata(1);
laps(1).start_idx = 1;
laps(1).direction = 0;
iLap = 1;

newUpThCross = 1;     % flag for new lap top search
newDownThCross = 1;     % flag for new lap top search
for i=1:length(Vdata)
    if Vdata(i) < th1    % search for min
        if newUpThCross       % start a new lap
            newUpThCross = 0;
            newDownThCross = 1;
            iLap = iLap + 1;
            laps(iLap).start_ts = Vts(i);
            laps(iLap).pos = Vdata(i);
            laps(iLap).start_idx = i;
            laps(iLap).direction = +1;
        end
        if Vdata(i) < laps(iLap).pos       % record new min if any
            laps(iLap).start_ts = Vts(i);
            laps(iLap).pos = Vdata(i);
            laps(iLap).start_idx = i;
        end
    end
    
    if Vdata(i) > th2   % search for max
        if newDownThCross       % start a new lap
            newUpThCross = 1;
            newDownThCross = 0;
            iLap = iLap + 1;
            laps(iLap).start_ts = Vts(i);
            laps(iLap).pos = Vdata(i);
            laps(iLap).start_idx = i;
            laps(iLap).direction = -1;
        end
        if Vdata(i) > laps(iLap).pos       % record new min if any
            laps(iLap).start_ts = Vts(i);
            laps(iLap).pos = Vdata(i);
            laps(iLap).start_idx = i;
        end
        
    end
    
end

% fix direction of first lap which was unknown above
% make first lap direction opposite of second lap's direction (laps alternate!)
laps(1).direction = -laps(2).direction;

% make sure laps cross the halfway point
middle=median(min(Vdata):max(Vdata));
i=1;
while true
    if any(laps(i).pos:laps(i+1).pos>middle) && any(laps(i).pos:laps(i+1).pos<middle)==0
        laps(i+1)=[];
    end
    i=i+1;
    if i>=length([laps.pos])
        if length([laps.pos])<iLap
            laps(1).direction = -laps(2).direction;
        end
        break
    end
end
end