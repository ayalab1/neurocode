function events = LoadEvents(filename)

%LoadEvents - Read events from file.
%
%  USAGE
%
%    events = LoadEvents(filename)
%
%    filename            event file name

% Copyright (C) 2004-2015 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

events.time = [];
events.description = [];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end

% Read file into cell array
file = fopen(filename,'r');
if file == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end
c = textscan(file,'%s','delimiter','\n');
fclose(file);

% Parse cell array (extract time and messages using regular expressions)
t = regexprep(c{1},'([^ \t]*).*','$1','once');
events.time = cellfun(@str2num,t);
events.description = regexprep(c{1},'[^ \t]*[ \t]*','','once');
%  t = regexprep(c{1},'([0-9]*[.]?[0-9]*[e]?[+-]?[0-9]*).*','$1','once');
%  events.description = regexprep(c{1},'[0-9]*[.]?[0-9]*[ \t]*','','once');


% Convert to seconds
if ~isempty(events.time), events.time(:,1) = events.time(:,1) / 1000; end
