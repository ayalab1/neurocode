function [histogram2d,h0,difference] = JointPETH(PETH1, PETH2,smooth)

% Example usage:
% [PETH1,t1] = PETH(ripples(:,2),deltas(:,2)); PETH1 = PETH1/mode(diff(t1)); % in Hz
% [PETH2,t2] = PETH(spindles(:,2),deltas(:,2)); PETH2 = PETH2/mode(diff(t2)); % in Hz
% [joint, expected, difference] = JointPETH(PETH1,PETH2,2);
% figure; PlotColorMap(joint,'x',t2,'y',t1); 
% xlabel('spindle rate centered on deltas'); ylabel('ripple rate centered on deltas'); 
% Note that the number of columns in PETH1 and PETH2 need to be the same (they are centered on the same events)


histogram2d = Smooth(PETH1'*PETH2,smooth);
h0 = Smooth(repmat(nanmean(PETH1),size(PETH1,1),1)'*repmat(nanmean(PETH2),size(PETH2,1),1),smooth);

histogram2d = histogram2d/size(PETH1,1); 
h0 = h0/size(PETH1,1);

histogram2d = sqrt(histogram2d); % make the final result in Hz, rather than Hz^2
h0 = sqrt(h0); % make the final result in Hz, rather than Hz^2
difference = histogram2d - h0;
