%% if you saved in the wrong format but you want to read digitalin 
%% step2 run
% D:\Documents\GitHub\neurocode\preProcessing\intan\Windows\IntanFileMerger.exe
% to merge all rhd 
  %step3
f = fopen('amplifier.dat','w'); % new file should be test.dat
fwrite(f,amplifier_data,'int16');
fclose(f);

% step4 in preprocessSession set fillMissingDatFiles true



%% deeal with digtal in
% step1 run read_Intan_RHD2000_file.m
% step2 directly extract digital_events_file
Nchan = size(board_dig_in_data,1); Nchan2 = Nchan +1;
clear pulses pulses2 digitalIn
for k = 1:Nchan 
    tester = board_dig_in_data(k,:);
    pulses{k} = strfind(tester,[0 1])';
    pulses2{k} = strfind(tester,[1 0])';
end
digital_on = pulses;
digital_off = pulses2;
disp('Done!');

for ii = 1:size(digital_on,2)
    if ~isempty(digital_on{ii})
        % take timestamp in seconds
        digitalIn.timestampsOn{ii} = digital_on{ii}/fs;
        digitalIn.timestampsOff{ii} = digital_off{ii}/fs;
        
        % intervals
        d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
        d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
        d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        digitalIn.ints{ii} = d;
        digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % duration
        
        clear intsPeriods
        intsPeriods(1,1) = d(1,1); % find stimulation intervals
        intPeaks =find(diff(d(1,:))>lag);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = d(2,intPeaks(jj));
            intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
        end
        intsPeriods(end,2) = d(2,end);  
        digitalIn.intsPeriods{ii} = intsPeriods;
    end
end

if exist('digitalIn')==1
    try save([sess.FileName '.DigitalIn.events.mat'],'digitalIn');
    catch
        save('digitalIn.events.mat','digitalIn');
    end
    keyboard
    
    clf
    for i=1:length(digitalIn.timestampsOn)
        intervals = digitalIn.intsPeriods{i};
        if ~isempty(intervals)
        PlotIntervals(intervals,'color','k','ylim',[0 1]+i-1,'alpha',1);
        end
    end
    ylim([0 length(digitalIn.timestampsOn)])
    mkdir('Pulses');
    saveas(gcf,'pulses\digitalIn.png')
else
    digitalIn = [];
end
%% step3 after preprecessSession
n=2;%the order of the session you are dealing with
j=1:size(board_dig_in_data,1),digitalIn{j} = [temp_digitalIn{j}+MergePoints.timestamps(n,1); digitalIn{j}];




