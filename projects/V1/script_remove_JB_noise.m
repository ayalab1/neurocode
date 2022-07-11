
rejectChannels = 1+[]; % for JB
nChannels = 160;
SSD_file = 'D:\test.dat';

basename = basenameFromBasepath(basepath);
datFile = [basepath,filesep, basename, '.dat'];
copyfile(datFile,SSD_file); 
disp([datestr(clock) ': dat file copied for ' basepath '. Finding noisy periods...']);
m = memmapfile(SSD_file, 'Format','int16','Writable',true);
data = reshape(m.data,nChannels,[]);
nSamples = size(data,2);
okChannels = ~ismember((1:size(data,1))',rejectChannels); okChannels = okChannels(1:4:end,:);
signal = mean(data(okChannels,:))';
bad = [false; abs(diff(signal))>200];
badIntervals = FindInterval(bad); badIntervals = [badIntervals(:,1)-1 badIntervals(:,2)+1];
badIntervals = ConsolidateIntervals(badIntervals,'epsilon',15);
save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
datestr((datenum(clock)))
noiseIntervalIndices = badIntervals;
noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;

disp([datestr(clock) ': Noise periods found. Removing noise in .dat file local copy for ' basepath '...']);
for i = 1:nChannels
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

disp([datestr(clock) ': Noise removed. Copying back local dat file to ' basepath '...']);
copyfile(SSD_file,datFile); 
clear m




















