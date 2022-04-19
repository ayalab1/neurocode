function [histogram2d,h0,difference] = JointPETH(PETH1, PETH2,smooth)

%JointPETH - produce a joint histogram for the co-occurrence of two sets of signals around events. 
% PETH1 and PETH2 should be in the format of PETH, see example usage below.
% This analysis tests for interactions. For example, the interaction of 
% ripples and spindles around the occurrence of delta waves. It is a good way
% to control whether the relationshops between two variables is entirely explained
% by a third variable (the events serving as basis for the PETHs). See Sirota et al. (2003)
%
% EXAMPLE
%
% [PETH1,t1] = PETH(ripples(:,2),deltas(:,2)); PETH1 = PETH1/mode(diff(t1)); % in Hz
% [PETH2,t2] = PETH(spindles(:,2),deltas(:,2)); PETH2 = PETH2/mode(diff(t2)); % in Hz
% [joint, expected, difference] = JointPETH(PETH1,PETH2,2);
% figure; PlotColorMap(joint,'x',t2,'y',t1); 
% xlabel('spindle rate centered on deltas'); ylabel('ripple rate centered on deltas'); 
% Note that the number of columns in PETH1 and PETH2 need to be the same (they are centered on the same events)
%
% Copyright (C) 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

histogram2d = Smooth(PETH1'*PETH2,smooth);
h0 = Smooth(repmat(nanmean(PETH1),size(PETH1,1),1)'*repmat(nanmean(PETH2),size(PETH2,1),1),smooth);

histogram2d = histogram2d/size(PETH1,1); 
h0 = h0/size(PETH1,1);

histogram2d = sqrt(histogram2d); % make the final result in Hz, rather than Hz^2
h0 = sqrt(h0); % make the final result in Hz, rather than Hz^2
difference = histogram2d - h0;
