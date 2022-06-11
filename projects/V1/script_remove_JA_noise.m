
rejectChannels = 1+[49 29 7 12 6 13 5 14 4 15 8 11 9 10]; % for Juan Antonio
nChannels = 96;

basename = basenameFromBasepath(basepath);
datFile = [basepath,filesep, basename, '.dat'];
m = memmapfile(datFile, 'Format','int16','Writable',true);
data = reshape(m.data,nChannels,[]);
nSamples = size(data,2);
okChannels = ~ismember((1:size(data,1))',rejectChannels); okChannels = okChannels(1:3:end,:);
signal = mean(data(okChannels,:))';
bad = [false; abs(diff(signal))>200];
badIntervals = FindInterval(bad); badIntervals = [badIntervals(:,1)-1 badIntervals(:,2)+1];
badIntervals = ConsolidateIntervalsFast(badIntervals,'epsilon',15);
save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
datestr((datenum(clock)))
noiseIntervalIndices = badIntervals;
noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;

for i = 1:nChannels
    i
    if any(rejectChannels==i),
        badIndices = sub2ind([nChannels,nSamples],i*ones(nSamples,1),(1:nSamples)');
        interpolated = 0;
    else
        badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
        goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
        badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
        goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
        goodValues = m.Data(goodIndices);
        interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
    end
    m.Data(badIndices) = int16(interpolated);
end



















