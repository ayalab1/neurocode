function createEVT(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used to create a .evt file from a generic ripple structure produced by
% various other analysis scripts. Note that this has also been adapted to
% pull timestamps from other non-ripple structures, so long as the naming
% convention (timestamps, peaks) is maintained. 
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
%   saveName      option for variability of the .evt file name. This will
%                 save as [basename].[saveName][#evt].evt (ie day6.R01.evt)
%   tVal          set as False to override the timestamp validation. This 
%                 is useful for when the user would like to input fewer 
%                 time arrays (ie no peak) and chooses the start time to 
%                 fill this slot. 
%   savePath      full path to save files to
%
% OUTPUT
%   None - message in command window will indicate successful formation of
%   the .evt file (with # of events for validation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize default values
p = inputParser;
addParameter(p,'basepath',pwd,@isstr); %basepath to access
addParameter(p,'unit',"s",@isstr); %unit for event epoch times
addParameter(p,'saveName', 'R', @isstr); %save name
addParameter(p,'tVal', 'True', @islogical); %use timestamp validation
addParameter(p,'savePath',pwd,@isstr); %place to save files

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
saveName = p.Results.saveName;
tVal = p.Results.tVal;
savePath = p.Results.savePath;

if (unit ~= "min")&&(unit ~= "s")&&(unit ~= "ms")
    error('Please enter a valid input for unit: [min, s, ms]');
end

if length(saveName) > 1
    error('saveName must be a single character');
end

% Assign start, end, and peak arrays accordingly if structure
if check_r
    swr_s = ripples.timestamps(:,1);
    swr_p = ripples.peaks;
    swr_e = ripples.timestamps(:,2);
end

% Validate inputs by checking timestamps are correctly ordered
if tVal
    if (swr_s(1)>=swr_p(1))||(swr_s(1)>=swr_e(1))||(swr_p(1)>=swr_e(1))
        error('Please validate that the timestamps are correctly ordered: start, peak, end');
    else
        disp('Input timestamps have been successfully validated');
    end
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
filename = basenameFromBasepath(basepath);

% Check if there is an existing .evt file in the current directory
oldPath = cd(savePath);
rippleFiles = dir(['*.' saveName '*.evt']);
cd(oldPath);
if isempty(rippleFiles)
    fileN = 1;
else
    % Set file index to next available value
    pat = ['.' saveName '[0-9].'];
    fileN = 0;
    for ii = 1:length(rippleFiles)
        token  = regexp(rippleFiles(ii).name,pat);
        val    = str2double(rippleFiles(ii).name(token+2:token+4));
        fileN  = max([fileN val]);
    end
    fileN = fileN + 1;
end

% Open the appropriately created new file
fid = fopen(sprintf(['%s%s%s.' saveName '%02d.evt'],savePath,filesep,filename,fileN),'w');

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