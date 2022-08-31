
function uLEDPulses = combineULEDPulses(varargin)
% combine analogPulses and digitalPulses an create a common structure
% containing uLEDPulses with its channel layout.
%
% INPUT:
%   basepath             Default pwd
%   analogPulses         Analog pulses with uLEd connection according to
%                           the following map
%   digitalPulses        Digital pulses with uLEd connection according to
%                           the following map
%   ledLayout            By default, see below
%   saveMat              Default, true.
%   force                Default, false;
%
% uLED map:
%  S1               S2               S3               S4
%  ________________________________________________________________________
%  L1: Analog Ch3   L1: Digit  Ch12  L1: Analog Ch6   L1: Digit  Ch15
%  L2: Digit  Ch11  L2: Analog Ch5   L2: Digit  Ch14  L2: Analog Ch8
%  L3: Analog Ch4   L3: Digit  Ch13  L3: Analog Ch7   L3: Digit  Ch16
%
% Legend:
%  uLED    Code
%__________________
%  S1L1     1
%  S1L2     2
%  S1L3     3
%  S2L1     4
%  S2L2     5
%  S2L3     6
%  S3L1     7
%  S3L2     8
%  S3L3     9
%  S4L1    10
%  S4L2    11
%  S4L3    12
%
% MV 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'analogPulses',[],@isstruct);
addParameter(p,'digitalPulses',[],@isstruct);
addParameter(p,'ledLayout',[],@iscell);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
analogPulses = p.Results.analogPulses;
digitalPulses = p.Results.digitalPulses;
ledLayout = p.Results.ledLayout;
saveMat = p.Results.saveMat;
force = p.Results.force;

prevPath = pwd;
cd(basepath);

targetFile = dir('*.uLEDPulses.event.mat');
if ~isempty(targetFile) && ~force
    disp('Pulses already sorted! Loading file.');
    load(targetFile.name);
    return
end


if isempty(ledLayout)
    ledLayout.channel =    [3 11  4 12  5 13  6 14  7 15  8 16];
    ledLayout.isAnalog =   [1  0  1  0  1  0  1  0  1  0  1  0];
    ledLayout.isDigital =  [0  1  0  1  0  1  0  1  0  1  0  1];
    ledLayout.code =       [1  2  3  4  5  6  7  8  9 10 11 12];
    ledLayout.shank =      [1  1  1  2  2  2  3  3  3  4  4  4];
    ledLayout.LED =        [1  2  3  1  2  3  1  2  3  1  2  3];
end


if isempty(analogPulses)
    analogPulses = getAnalogPulses('analogCh',ledLayout.channel(ledLayout.isAnalog==1));
end

if isempty(digitalPulses)
    digitalPulses = getDigitalIn;
end

%% Collect pulses
timestamps = []; code = []; shank = []; LED = []; pulsesNumber = [];
for ii = 1:length(ledLayout.channel)
        if ledLayout.isAnalog(ii)
            timestamps = [timestamps; analogPulses.timestamps(analogPulses.analogChannel==ledLayout.channel(ii),:)];
            code =       [code      ; ledLayout.code(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            shank =      [shank     ; ledLayout.shank(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            LED =        [LED       ; ledLayout.LED(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            pulsesNumber(ii) = length(find(analogPulses.analogChannel==ledLayout.channel(ii)));
        elseif ledLayout.isDigital(ii)
            timestamps = [timestamps; digitalPulses.ints{ledLayout.channel(ii)}'];
            code =       [code      ; ledLayout.code(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},2),1)];
            shank =      [shank     ; ledLayout.shank(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},2),1)];
            LED =        [LED       ; ledLayout.LED(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},2),1)];
            pulsesNumber(ii) = size(digitalPulses.ints{ledLayout.channel(ii)},2);
        end
end

[~, idx] = sort(timestamps(:,1));
timestamps = timestamps(idx,:);
code = code(idx);
shank = shank(idx);
LED = LED(idx);


%% Generate output
uLEDPulses.timestamps = timestamps;
uLEDPulses.code = code;
uLEDPulses.shank = shank;
uLEDPulses.LED = LED;
uLEDPulses.pulsesNumber = pulsesNumber;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.uLEDPulses.events.mat'],'uLEDPulses');
end

cd(prevPath);
end