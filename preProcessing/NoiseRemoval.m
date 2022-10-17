function NoiseRemoval(basepath)

%[NoiseRemoval] - [removes noise periods in LFP/DAT files]
%
%  [removes noisy periods in LFP/DAT files, ensuring that noisy periods are
%  fully removed]
%
%
%  INPUTS
%
%    [basepath]      [basepath of session with DATs to have noise removed. 
%                    Note that basepath must be provided]
%
%
%  OUTPUTS
%    
%    NA
%
%  NOTE
% 
%   [uses helper function - 
%       function y = linspaceVector(d1, d2)
%       function Logical = Unfind(indices, n)
%       function s = CumSum(data,stops)]
%    HLR - are ^ needed since these should be in our "utilities" subfolder?
%
%
%   TODO: added options to make basepath "pwd"
%
%  SEE ALSO
%
% [KM and 2022 Ralitsa Todorova] [2021-2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%------------------------------------------------------------------------



session = getSession('basepath',basepath);
datFile = [basepath,filesep, session.general.name, '.dat'];
datSamplingRate = session.extracellular.sr;

warning(['This function will change the file ' datFile '.']);
nChannels = session.extracellular.nChannels;
lfpStructure = getLFP('all','basepath',basepath);
lfpSamplingRate = lfpStructure.samplingRate;

% Identify Noise intervals
lfp = mean(lfpStructure.data,2);

[b,a] = cheby2(5,20,300/(lfpSamplingRate/2),'low');
slow = filtfilt(b,a,lfp);
filtered = lfp-slow; % this is the signal filtered in the 300Hz+ frequency band
power = abs(hilbert(filtered));

% Find noise outliers:
q1 = quantile(power,0.25);
q3 = quantile(power,0.75);
d = q3-q1;
threshold = q3+3*d; % this Tukey's definition of outliers that are "far out" (normal outliers would be at "q3+1.5*d" but we want to be conservative here)
noiseIntervals = FindInterval(power>threshold);
% add at least 2.5 ms before and after each interval
noiseIntervals = bsxfun(@plus,noiseIntervals,[-2.5 +2.5]/1000*lfpSamplingRate);

% the units of "noiseIntervals" are in lfp indices. Convert them to dat indices
noiseIntervalIndices = [floor(noiseIntervals(:,1)/lfpSamplingRate*datSamplingRate) ceil(noiseIntervals(:,2)/lfpSamplingRate*datSamplingRate)];

% Merge noise intervals within "epsilon" of each other. This step is important 
% (even with epsilon=0s) because it ensures the points immediately preceding 
% and immediately following the boundaries (e.g. noiseIntervalIndices(j,1)-1 
% and noiseIntervalIndices(j,2)+1) are noise-free and we can interpolate from them
noiseIntervalIndices = ConsolidateIntervals(noiseIntervalIndices,'epsilon',5/1000*datSamplingRate); % epsilon = 5ms

warning(['Removing a total of ' num2str(sum(diff(noiseIntervalIndices,[],2))/datSamplingRate) 's of noisy data from the .dat file']);
%%
m = memmapfile(datFile, 'Format','int16','Writable',true);
nSamples = round(length(m.Data)/nChannels);
% make sure that the very first and the very last points in the dat-file are not selected to be removed (this would make interpolation impossible)
noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;

try
    noiseIntervals = noiseIntervalIndices/datSamplingRate;
    disp(['Saving noise intervals in ' fullfile(basepath,'noiseIntervals.events.mat')]);
    save(fullfile(basepath,'noiseIntervals.events.mat'),'noiseIntervals');
end
for i = 1:nChannels
    badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
    goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
    badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
    goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
    goodValues = m.Data(goodIndices);
    interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
    m.Data(badIndices) = int16(interpolated);
end


% ------------------------------- Helper functions -------------------------------

function y = linspaceVector(d1, d2)

% It's like linspace, but d1 and d2 and vectors
% It's a faster equivalent to a loop calling linscape for each [d1 d2] pair:
% for i=1:length(d1),y = [y;linspace(d1(i),d2(i))']; end
%
% Example:
% linspaceVector([2;5;23],[2;8;23]) = [2;5;6;7;8;23]
%
% Copyright (C) 2020 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

n = d2-d1;

% repeat original value n times
y0 = repelem(d1,n+1);
nRows = size(y0,1);

% add 1 for each additional value after d1
indicesOfOriginalValues = cumsum([0;n(1:end-1)]+1);
originalValues = Unfind(indicesOfOriginalValues,nRows);
toAdd = CumSum(ones(nRows,1),originalValues)-1;

y = y0 + toAdd;


function Logical = Unfind(indices, n)

% When you have a logical vector, 'find' is useful as it gives you the
% indices of the non-zero values.
% This function performs the opposite operation, giving you a logical
% vector with ones in the positions of the indices provided.
%
% If you want the length of the resulting logical to be different from the
% last index provided, provide n (default = last index).
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
 
if nargin<2,
    n=max(indices);
end
 
Logical = false(n,1);
Logical(indices) = true;

function s = CumSum(data,stops)

%CumSum - Cumulative sum of elements. Partial sums can also be computed.
%
%  USAGE
%
%    sum = CumSum(data,stops)
%
%    data           data to sum
%    stops          optional logical indices where sum should be restarted
%
%  EXAMPLE
%
%    % Simple cumulative sum
%    s = CumSum([1 4 6 2 13]);
%
%    % Partial cumulative sums
%    s = CumSum([1 4 6 2 13 2 4 6 5 10 1],[1 0 0 0 0 1 0 0 0 0 0])
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


data = data(:);
if nargin == 2,
	stops = logical(stops(:));
end

% Simple cumulative sum
s = cumsum(data);
if nargin == 1, return; end

% Use stops to restart cumulative sum (tricky vector computation)
stops(1) = 0;
i = find(stops);
k = s(i-1);
dk = diff([0;k]);
z = zeros(size(data));
z(i) = dk;
s = s-cumsum(z);


