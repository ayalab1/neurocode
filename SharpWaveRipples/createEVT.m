function createEVT(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used to create a .evt file from a generic ripple structure produced by
% various other analysis scripts
% 
% 
% INPUT
%   ripples       structure which must include the labeled timestamps of
%                 the start, stop, and peak power of each ripple event
%
%       OR
%   
%   start         array of start times for each ripple event
%   peak          array of times denoting peak power of each ripple event
%   stop          array of end times for each ripple event
%
%
% OPTIONS
%   basepath      path to a single session to run
%   unit          unit of time which is used for these timestamps. Valid 
%                 inputs include [min, s, or ms]
%
% OUTPUT
%   None - message in command window will indicate successful formation of
%   the .evt file (with # of events for validation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize default values
p = inputParser;
addParameter(p,'basepath',pwd,@isstr); %basepath to access
addParameter(p,'unit',"s",@isstr); %unit for event epoch times

% Determine input type and require inputs accordingly
if isstruct(varargin{1})
    addRequired(p, 'ripples', @isstruct);
    parse(p,varargin{:});
    ripples = p.Results.ripples;
    check_r = 1;
else
    addRequired(p, 'start', @isnumeric);
    addRequired(p, 'peak', @isnumeric);
    addRequired(p, 'stop', @isnumeric);
    parse(p, varargin{:});
    swr_s = p.Results.start;
    swr_p = p.Results.peak;
    swr_e = p.Results.stop;
    check_r = 0;
end

% Assign parameters to variables accordingly
basepath = p.Results.basepath;
unit = p.Results.unit;
if (unit ~= "min")&&(unit ~= "s")&&(unit ~= "ms")
    error('Please enter a valid input for unit: [min, s, ms]');
end

% Assign start, end, and peak arrays accordingly if structure
if check_r
    swr_s = ripples.timestamps(:,1);
    swr_p = ripples.peaks;
    swr_e = ripples.timestamps(:,2);
end

% .evt file should be in ms - convert accordingly
if unit == "s"
    swr_s = swr_s.*1000;
    swr_p = swr_p.*1000;
    swr_e = swr_e.*1000;
elseif unit == "min"
    swr_s = swr_s.*(60*1000);
    swr_p = swr_p.*(60*1000);
    swr_e = swr_e.*(60*1000);
end

% Set file name standards
Filebase = basenameFromBasepath(basepath);
Filebase = fullfile(basepath,Filebase);
[pathname, filename, extname] = fileparts(Filebase);
if isempty(pathname)
    pathname = pwd;
end

% Check if there is an existing .evt file in the current directory
rippleFiles = dir('*.R*.evt');
if isempty(rippleFiles)
    fileN = 1;
else
    % Set file index to next available value
    pat = '.R[0-9].';
    fileN = 0;
    for ii = 1:length(rippleFiles)
        token  = regexp(rippleFiles(ii).name,pat);
        val    = str2double(rippleFiles(ii).name(token+2:token+4));
        fileN  = max([fileN val]);
    end
    fileN = fileN + 1;
end

% Open the appropriately created new file
fid = fopen(sprintf('%s%s%s.R%02d.evt',pathname,filesep,filename,fileN),'w');

fprintf(1,'Writing event file ...\n');

% Write to the file
for i = 1:size(swr_s, 1)
    fprintf(fid,'%9.1f\tstart\n',swr_s(i));
    fprintf(fid,'%9.1f\tpeak\n',swr_p(i));
    fprintf(fid,'%9.1f\tstop\n',swr_e(i));
end
fclose(fid);
disp(['Event file successfully created - contains ' num2str(size(swr_s,1)) ' events.']);

end