function [joint,expected,difference] = JointPETH(PETH1, PETH2,smooth)

%JointPETH - produce a joint histogram for the co-occurrence of two sets of signals around events. 
% PETH1 and PETH2 should be in the format of PETH, see example usage below.
% This analysis tests for interactions. For example, the interaction of 
% ripples and spindles around the occurrence of delta waves. It is a good way
% to control whether the relationships between two variables is entirely explained
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
% Note: sometimes the difference between "joint" and "expected" may be dominated due to 
% brain state effects (e.g. if both ripples are spindles are more common around delta
% waves taking place in early SWS and have decreased rates around delta waves in late
% SWS, then all the values of "joint" would be larger than the value of "expected". 
% In such a case, to investigate the timing effects in particular and ignore such
% global changes (correlations across the rows of "PETH1" and "PETH2"), consider 
% normalizing the rows of the PETHs before calling JointPETH.
% e.g. nPETH1 = PETH1./sum(PETH1,2); nPETH2 = PETH2./sum(PETH2,2);
% [joint,expected,difference] = JointPETH(nPETH1, nPETH2,smooth)
%
% Copyright (C) 2018-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

joint = Smooth(PETH1'*PETH2,smooth);
expected = Smooth(repmat(nanmean(PETH1),size(PETH1,1),1)'*repmat(nanmean(PETH2),size(PETH2,1),1),smooth);

joint = joint/size(PETH1,1); 
expected = expected/size(PETH1,1);

joint = sqrt(joint); % make the final result in Hz, rather than Hz^2
expected = sqrt(expected); % make the final result in Hz, rather than Hz^2
difference = joint - expected;
