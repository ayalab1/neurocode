function ProcessBinary(inputName,outputName,nChannels,f,varargin)

%ProcessBinary - Process binary data file.
%
% Process binary data file, e.g. filter LFP file. This function loads the
% binary file segment by segment, calls a user-supplied function to process
% each data segment, and optionally saves the result to a new file.
%
% If the function requires overlapping segments (e.g. to avoid edge effects),
% it will receive data in the form [o1;s;o2], where s is the segment to process,
% and o1 and o2 are portions of the previous and next segments (half the overlap
% size each). The function should return the processed segment, WITHOUT the
% overlapping portions.
%
%  USAGE
%
%    ProcessBinary(inputName,outputName,nChannels,f,<options>)
%
%    inputName      binary input file
%    outputName     binary output file (optional, see below)
%    nChannels      number of channels in the input file
%    f              function handle
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'parameters'  additional parameters for f (cell array)
%     'overlap'     overlap in # samples (default = 0)
%     'segment'     segment length in # samples (default = 2^16)
%    =========================================================================
%
%  EXAMPLES
%
%    % Change sign (using an anonymous function)
%    ProcessBinary('input.dat','output.dat',1,@(x) -x);
%
%    % Low-pass filter and square using the following funtion:
%    %   function y = CustomFilter(x,b,a)
%    %     y = filtfilt(b,a,x).^2;
%    %     y = y(251:end-250);
%    ProcessBinary('input.dat','output.dat',1,@CustomFilter,'parameters',{b,a},'overlap',500);
%
%    % If you require more elaborate functionality than just saving the processed
%    % segment (e.g. save several files, or use a custom file format), pass an
%    % empty output filename and include the appropriate code in your function.
%    % For instance,
%    %   function y = SaveMax(x,f)
%    %     m = max(x);
%    %     fwrite(f,m);
%    ProcessBinary('input.dat','',1,@SaveMax,'parameters',{'output.dat'});
%    

% Copyright (C) 2004-2018 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
nOverlap = 0;
segmentLength = 2^16;
parameters = {};

% Check parameters
if nargin < 4,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'overlap',
			nOverlap = varargin{i+1};
			if ~isiscalar(nOverlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
			end
			if rem(nOverlap,2) ~= 0,
				error('Overlap must be even (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
			end
		case 'segment',
			segmentLength = varargin{i+1};
			if ~isiscalar(segmentLength,'>0'),
				error('Incorrect value for property ''segment'' (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
			end
		case 'parameters',
			parameters = varargin{i+1};
			if ~iscell(parameters),
				error('Incorrect value for property ''parameters'' (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).']);
	end
end

% Open input and output files
inputFile = fopen(inputName,'r');
if ~isempty(outputName), outputFile = fopen(outputName,'w'); end

% Process first segment
if nOverlap == 0,
	overlap = zeros(nChannels,0);
else
	% Flip first segment LR and prepend it to avoid edge effects
	overlap = fread(inputFile,[nChannels,nOverlap/2],'int16');
	overlap = fliplr(overlap);
	frewind(inputFile);
end
segment = fread(inputFile,[nChannels,segmentLength],'int16');
segment = [overlap,segment]';
processed = feval(f,segment,parameters{:});
if ~isempty(outputName), fwrite(outputFile,processed,'int16'); end
overlap = segment(end-(nOverlap-1):end,:);

% Process subsequent segments
while ~feof(inputFile),
	segment = fread(inputFile,[nChannels,segmentLength],'int16');
	segment = [overlap;segment'];
	processed = feval(f,segment,parameters{:});
	if ~isempty(outputName), fwrite(outputFile,processed,'int16'); end
	overlap = segment(end-(nOverlap-1):end,:);
end

% Process trailing unprocessed segment
if nOverlap ~= 0,
	% Flip last segment UD and append it to avoid edge effects
	segment = segment(end-(nOverlap-1):end,:);
	tail = flipud(segment(end-(nOverlap/2-1):end,:));
	segment = [segment;tail];
	processed = feval(f,segment,parameters{:});
	if ~isempty(outputName), fwrite(outputFile,processed,'int16'); end
end

% Close input and output files
fclose(inputFile);
if ~isempty(outputName), fclose(outputFile); end

