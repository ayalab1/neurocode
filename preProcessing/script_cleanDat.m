
rejectChannels = 1+[]; 
SSD_file = 'D:\current_cleaning.dat';

basename = basenameFromBasepath(basepath);
session = getStruct(basepath,'session');
nChannels = sum(cellfun(@numel,session.extracellular.spikeGroups.channels));

timeAroundArtefactToRemove = 0; % in seconds
mergeArtefactsCloserThan = 5; % in seconds
smooth = 20000;
threshold = 30;
peakThreshold = 35;

datFile = [basepath,filesep, basename, '.dat'];
copyfile(datFile,SSD_file);
disp([datestr(clock) ': dat file copied for ' basepath '. Finding noisy periods...']);
m = memmapfile(SSD_file, 'Format','int16','Writable',true);
data = reshape(m.data,nChannels,[]);
nSamples = size(data,2);
% Find noisy periods
if exist(fullfile(basepath,'noiseIntervalsDat.events.mat'),'file')
    load(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
else
    okChannels = ~ismember((1:size(data,1))',rejectChannels); okChannels = okChannels(1:4:end,:);
    signal = mean(data(okChannels,:))';
    d = [0;abs(diff(signal))]; smoothed = Smooth(d,smooth);
    bad = smoothed>threshold;
    badIntervals = FindInterval(bad); 
    badIntervals = ConsolidateIntervalsFast(badIntervals,'epsilon',ceil(mergeArtefactsCloserThan*20000));
    badIntervals(diff(badIntervals,[],2)==0,:) = []; % ignore single bin intervals
    badIntervals = [badIntervals(:,1)-ceil(timeAroundArtefactToRemove*20000) badIntervals(:,2)+ceil(timeAroundArtefactToRemove*20000)];
    badIntervals = ConsolidateIntervalsFast(badIntervals);
    badIntervals(badIntervals<2) = 2; badIntervals(badIntervals>nSamples-1) = nSamples-1;
    peaks = zeros(size(badIntervals(:,1)));
    for i=1:size(badIntervals,1)
        peaks(i) = max(smoothed(badIntervals(i,1):badIntervals(i,2)));
    end
    badIntervals(peaks<peakThreshold,:) = [];

    t = (1:length(smoothed))/20000;
    figure; plot(t(1:1000:end)',smoothed(1:1000:end));
    PlotIntervals(t(badIntervals),'color','k');
    drawnow
    choice = str2double(input('Are you happy with this detection of noisy intervals? (1=yes, DELETE THIS FROM MY .DAT FILE; 0=no, keep groups as they are; -1=enter debug mode): ','s'));

    if choice==0
        return
    elseif choice == -1
        keyboard
    end

    save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
end
datestr((datenum(clock)))
noiseIntervalIndices = badIntervals;
noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;

disp([datestr(clock) ': Noise periods found. Removing noise in .dat file local copy for ' basepath '...']);
for i = 1:nChannels
    badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
    goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
    badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
    goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
    goodValues = m.Data(goodIndices);
    interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
    m.Data(badIndices) = int16(interpolated);
end

disp([datestr(clock) ': Noise removed. Copying back local dat file to ' basepath '...']);
copyfile(SSD_file,datFile); 
clear m





















