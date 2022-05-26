
rejectChannels = []; % for AO52
rejectChannels = 1+[8 12 14 15]; % for AO52
nChannels = 64;

basename = basenameFromBasepath(basepath);
datFile = [basepath,filesep, basename, '.dat'];
m = memmapfile(datFile, 'Format','int16','Writable',true);
data = reshape(m.data,64,[]);
nSamples = size(data,2);
okChannels = ~ismember((1:size(data,1))',rejectChannels);
signal = mean(data(okChannels,:))';
bad = [false; abs(diff(signal))>200];
badIntervals = FindInterval(bad); badIntervals = [badIntervals(:,1)-1 badIntervals(:,2)+1];
badIntervals = ConsolidateIntervals(badIntervals,'epsilon',15);
save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
datestr((datenum(clock)))
noiseIntervalIndices = badIntervals;
noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;
for i = 1:nChannels
    i
    badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
    goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
    badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
    goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
    goodValues = m.Data(goodIndices);
    interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
    m.Data(badIndices) = int16(interpolated);
end

%% For Jean (doesn't have an effect on cluster quality, do not apply)
% 
% rejectChannels = 1+[22 13 15 47 19 51 49 57 55 43 54 39 50 52 41 48 37]; % for Jean
% rejectChannels = 1+[22 13 15 24 47 19 51 49 57 55 43 54 39 50 52 41 48 37]; % for Jean\day54
% nChannels = 64;
% 
% basename = basenameFromBasepath(basepath);
% datFile = [basepath,filesep, basename, '.dat'];
% m = memmapfile(datFile, 'Format','int16','Writable',true);
% data = reshape(m.data,64,[]);
% nSamples = size(data,2);
% okChannels = ~ismember((1:size(data,1))',rejectChannels);
% signal = mean(data(okChannels,:))';
% bad = [false; abs(diff(signal))>200];
% d = [0;diff(signal)];
% previous = d(1:end-1); current = d(2:end);
% up = previous < 0 & current > 0;up(end+1) = 0;
% down = previous > 0 & current < 0;down(end+1) = 0;
% downIntervals = find(down); downIntervals = [downIntervals(1:end-1) downIntervals(2:end)];
% upIntervals = find(up); upIntervals = [upIntervals(1:end-1) upIntervals(2:end)];
% indicesBad = find(bad);
% indicesBadDownIntervals = unique(FindClosest(downIntervals(:,1),indicesBad,'lower'));
% indicesBadUpIntervals = (FindClosest(upIntervals(:,1),indicesBad,'lower'));
% badIntervals = sortrows([downIntervals(indicesBadDownIntervals,:); upIntervals(indicesBadUpIntervals,:)]);
% badIntervals = ConsolidateIntervals(badIntervals,'epsilon',15);
% badIntervals = (badIntervals(diff(badIntervals,[],2)<400,:)); % take only the ones shorter than 20ms
% save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
% datestr((datenum(clock)))
% noiseIntervalIndices = badIntervals;
% noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;
% for i = 1:nChannels
%     i
%     badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
%     goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
%     badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
%     goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
%     goodValues = m.Data(goodIndices);
%     interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
%     m.Data(badIndices) = int16(interpolated);
% end
























