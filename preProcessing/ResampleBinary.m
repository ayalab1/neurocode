function ResampleBinary(inputName,outputName,nChannels,up,down)

%ResampleBinary - Resample binary data file.
%
% Resample binary data file, e.g. create LFP file from raw data file.
%
%  USAGE
%
%    ResampleBinary(inputName,outputName,nChannels,up,down)
%
%    inputName      binary input file
%    outputName     binary output file
%    nChannels      number of channels in the file
%    up             upsampling integer factor
%    down           downsampling integer factor
%
%  NOTE 1
%
%    This function is provided for convenience. It simply calls <a href="matlab:help ProcessBinary">ProcessBinary</a>
%    using the same parameters. See this function for details.
%
%  NOTE 2
%
%    The actual resampling ratio is up/down.
%
%    Here is a list of typical values for Spike2 recording systems:
%
%    FROM           TO       UP     DOWN
%    =====================================
%    20000          1025     1      16
%    19531.25       20000    128    125
%    19531.25       1025     128    125*16
%    19841.29...    20000    126    125
%    19841.29...    1025     126    125*16
%    20284          20000    5000   5071
%    20284          1025     5000   5071*16

% Copyright (C) 2004-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin ~= 5,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ResampleBinary">ResampleBinary</a>'' for details).');
end

segmentLength = 2^16  - mod(2^16,down);
% Number of overlapping points per channel, chosen so that both resampled and original overlaps are integers
resampledOverlap = 8*up;
originalOverlap = resampledOverlap * down/up;

ProcessBinary(inputName,outputName,nChannels,@ResampleSegment,'parameters',{up,down},'overlap',originalOverlap,'segment',segmentLength);

function y = ResampleSegment(x,up,down)

resampledOverlap = 8*up;
y = resample(x,up,down);
y = y(resampledOverlap/2+1:end-resampledOverlap/2,:)';
