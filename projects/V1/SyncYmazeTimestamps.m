% To synchronize matlab-saved matrix and digital inputs, we'll assume that
% both clocks are correct and look for a particular delay between the two:

a = dlmread('arduino-outputs-20220408T134338.times');
matrix = dlmread('task-ymaze-20220408T134338.times');
load('digitalIn.events.mat','digitalInp');

aa = []; for i=1:3, qq = a(strfind(a(:,i+1)',[0 1])+1)'; qq(:,2) = i; aa = [aa; qq]; end

q = digitalInp.timestampsOn;
g = Group(q{[1 3 4]});

delay0 = aa(1)-g(1);
[h,ht] = mPETH(aa(:,1),g(:,1)+delay0,'durations',[-1 1]*100,'nBins',10001);
delay1 = delay0 + ht(findmax(h));
[h,ht] = mPETH(aa(:,1),g(:,1)+delay1,'durations',[-1 1]*1,'nBins',1001);
delay2 = delay1 + ht(findmax(h));

% this is the delay attempting to place the digital inputs and the first MATLAB detection as simultanous. 
% In reality, matlab samples the arduino every ~80ms. So the PETH of matlab-recorded signals should be
% high from 0 to 80ms (depending on when in the window the signal was sampled).
% Find the perfect 80ms window around the digital inputs, maximizing the coincidences:
possibleDelays = ((-80:1:80))'/1000; possibleDelays(:,2) = possibleDelays+80/1000;
count = nan(size(possibleDelays,1),1);
for i=1:size(possibleDelays,1)
    count(i,1) = sum(CountInIntervals(aa(:,1),g(:,1)+possibleDelays(i,:) + delay2));
end
delay = possibleDelays(findmax(count)) + delay2;

ymaze = matrix;
ymaze(1) = 0;
ymaze(:,1) = ymaze(:,1)-delay;

% dlmwrite('ymaze_synced.times',ymaze,'precision','%.6f');