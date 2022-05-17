
function [pulses] = getAnalogPulses(varargin)
% [pul, val, dur] = getAnalogPulses(varargin)
%
% Find square pulses. If not argument it provide, it tries to find pulses
% in intan analog-in file.
%
% <OPTIONALS>
% analogCh      List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
% data          R x C matrix with analog data. C is data, R should be
%               greater than 1.
% samplingRate            Sampling frequency (in Hz), default 20000.
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 20s)    
% filename      File to get pulses from. Default, data file with folder
%               name in current directory
% manualThr     Check manually threslhold amplitude (default, false)
% groupPulses   Group manually train of pulses (default, false)
% basepath      Path with analog data files to get pulses from.
% minDur        pulses with shorter duration than this are removed
%
% OUTPUTS
%               pulses - events struct with the following fields
% timestamps    C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses. 
% amplitude     values of the pulses with respect baleline (normalized as 0).
% duration      Duration of the pulses. Note that default samplingRate is 20000.
% eventID       Numeric ID for classifying various event types (C X 1)
% eventIDlabels label for classifying various event types defined in eventID (cell array C X 1)  
% intsPeriods   Stimulation periods, as defined by perioLag
%
% Manu-BuzsakiLab 2018
% Antonio FR, 10/21

%% Parse options
p = inputParser;
addParameter(p,'analogCh',[],@isnumeric);
addParameter(p,'data',[],@isnumeric);
addParameter(p,'samplingRate',20000,@isnumeric);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'filename',[],@isstring);
addParameter(p,'periodLag',20,@isnumeric);
addParameter(p,'manualThr',false,@islogical);
addParameter(p,'groupPulses',false,@islogical);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'useGPU',true,@islogical);
addParameter(p,'minDur',[],@isnumeric);

parse(p, varargin{:});
samplingRate = p.Results.samplingRate;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
manualThr = p.Results.manualThr;
d = p.Results.data;
analogCh = p.Results.analogCh;
groupPulses = p.Results.groupPulses;
basepath = p.Results.basepath;
useGPU = p.Results.useGPU;
minDur = p.Results.minDur;

prevPath = pwd;
cd(basepath);

%%
filetarget = split(pwd,filesep); filetarget = filetarget{end};
if exist([filetarget '.pulses.events.mat'],'file') 
    disp('Pulses already detected! Loading file.');
    load([filetarget '.pulses.events.mat']);
    if ~isempty(analogCh) && isnumeric(analogCh)
        maskPulses = ismember(pulses.analogChannel, analogCh);
        pulses.timestamps = pulses.timestamps(maskPulses,:);
        pulses.amplitude = pulses.amplitude(maskPulses,:);
        pulses.duration = pulses.duration(maskPulses,:);
        pulses.eventGroupID = pulses.eventGroupID(maskPulses,:);
        pulses.analogChannel = pulses.analogChannel(maskPulses,:);
        
    end
    return
end

analogFile = []; IntanBuzEd = [];
f=dir('analogin.dat');                                                     % check file
if isempty(f) || f.bytes == 0                                              % if analogin is empty or doesn't exist
    warning('analogin.dat file is empty or does not exist, was the recording made in Intan Buzsaki edition?');
    
    f = dir('*amplifier*.dat');                                        % is exist amplifier
    if isempty(f)
        analogFile = split(pwd,'\'); analogFile = analogFile{end};
        analogFile = strcat(analogFile,'.dat');
    else
        analogFile = f.name;
    end
    
    if isempty(analogCh)
        error('No posible to run getAnalogPulses from Intan Buzsaki Ed with no analogCh inputs!');
    else 
        analogCh = analogCh+1; % 0 to 1 index
    end
    parameters = LoadParameters(pwd); % read xml
    samplingRate = parameters.rates.wideband;
    nChannels = parameters.nChannels;
    IntanBuzEd = 1;
else
    disp('Using analogin.dat...');
    analogFile = 'analogin.dat';  
    try [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
        board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
        read_Intan_RHD2000_file_bz;
    catch
        disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
    end
    if ~exist('board_adc_channels')
        disp('Trying intan metadata from recording folder...');
        filetarget = split(filetarget,'_'); filetarget = filetarget{1};
        localPaths = dir([filetarget, '*']);
        cd(localPaths(1).name);
        [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
        board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
        read_Intan_RHD2000_file_bz;
        cd ..
    end
    samplingRate = frequency_parameters.board_adc_sample_rate;
    nChannels = length(board_adc_channels); % ADC input info from header file
    if isempty(analogCh)
        analogCh = 1:nChannels;
    end
    fileTargetAnalogIn =  dir('analogin*.dat');
    mAnalogIn = memmapfile(fullfile(basepath,fileTargetAnalogIn.name),'Format','uint16','Writable', true);
    dataAnalogIn = reshape(mAnalogIn.Data,size(board_adc_channels,2),[]);
    IntanBuzEd = 0;
end

h=figure;
% set(gcf,'Position',[100 -100 2500 1200]);
for jj = 1 : length(analogCh)
    fprintf(' ** Channel %3.i of %3.i... \n',jj, length(analogCh));
    disp('    Loading file...');
    if IntanBuzEd
        d = bz_LoadBinary(analogFile, 'frequency', samplingRate, 'nChannels', nChannels,'channels', analogCh(jj));
    else
        d = dataAnalogIn(analogCh(jj),:);
    end    
    xt = linspace(1,length(d)/samplingRate,length(d));
    
    if any(d<0) % if signal go negative, rectify
        d = d - min(d);
    end
    
    if ~manualThr
        thr = 150*median(d(1:100:end)/0.6745); % computing threshold
        if thr == 0 || ~any(d>thr)
            disp('Trying 5*std threshold...');
            thr = 4.5 * std(double(d));
        end
    else
        h2 = figure;
        plot(xt(1:100:end), d(1:100:end));
        xlabel('s'); ylabel('amp');
        title('Select threshold with the mouse and press left click...');
        [~,thr] = ginput(1);
        hold on
        plot([xt(1) xt(end)],[thr thr],'-r');
        pause(1);
        close(h2);
    end
    
    eventGroup = [];
    if groupPulses
        h = figure;
        plot(xt(1:100:end), d(1:100:end));
        hold on
        xlabel('s'); ylabel('amp');
        title('Group stimulation periods by pressing left click. Press enter when done.');
        selecting = 1;
        while selecting
            [x,~] = ginput(2);
            if ~isempty(x)
                plot([x(1) x(2)],[thr thr]);
                eventGroup = [eventGroup; x'];
            else
                selecting = 0;
            end
        end        
    end
    
    disp('    Thresholding signal and finding pulses...');
    dBin = (d>thr); % binarize signal
    locsA = find(diff(dBin)==1)/samplingRate; % start of pulses
    locsB = find(diff(dBin)==-1)/samplingRate; % end of pulses
    
    if useGPU
        temp = gpuArray(locsA(1:min([length(locsA) length(locsB)])));
    else
        temp = (locsA(1:min([length(locsA) length(locsB)])));
    end
    if size(temp,1) > size(temp,2)
        temp = temp';
    end

    temp2 = temp;
    temp(2,:) = 0;
    parfor ii = 1 : length(temp2) % pair begining and end of the pulse
        try temp(2,ii) =  locsB(find(locsB - temp2(ii) ==...
            min(locsB(locsB > temp2(ii)) - temp2(ii))));
        catch
            keyboard;
        end
    end
    temp(:,find(temp(1,:) == 0)) = [];
    
    pul{jj} = gather(temp);
    clear temp temp2
    
    baseline_d = int32(median(d(1:100:end)));
    val{jj}=[];
    for ii = 1 : size(pul{jj},2) % value of the pulse respect basaline
        val{jj}(ii) = median(int32(d(int32(pul{jj}(1,ii) * samplingRate : pul{jj}(2,ii) * samplingRate)))) - baseline_d;
    end
    
    pul{jj} = pul{jj} - offset;
    % discard pulses < 2 * median(abs(x)/0.6745) as noise or pulses in negatives times
    idx = find((val{jj} < thr*0.4) | pul{jj}(1,:)<0);
    val{jj}(idx) = [];
    pul{jj}(:,idx) = [];
    
    if ~isempty(pul{jj})
        dur{jj} = pul{jj}(2,:) - pul{jj}(1,:); % durantion
        
        stimPer{jj}(1,1) = pul{jj}(1,1); % find stimulation intervals
        intPeaks =find(diff(pul{jj}(1,:))>lag);
        for ii = 1:length(intPeaks)
            stimPer{jj}(ii,2) = pul{jj}(2,intPeaks(ii));
            stimPer{jj}(ii+1,1) = pul{jj}(1,intPeaks(ii)+1);
        end
        stimPer{jj}(end,2) = pul{jj}(2,end);
    else
        dur{jj} = [];
        stimPer{jj} = [];
    end
    
    eventGroupID{jj} = ones(size(dur{jj})) * jj;
    if ~isempty(eventGroup)
        for kk = 1:size(eventGroup,1)
            eventGroupID{jj}(pul{jj}(1,:) >= eventGroup(kk,1) & pul{jj}(1,:) <= eventGroup(kk,2)) = jj + size(d,1) + kk - 2;
        end
    end
    eventChannel{jj} = ones(size(dur{jj})) * analogCh(jj);
    
    % d = gpuArray(d); locsA = gpuArray(locsA);
    figure(h);
    subplot(length(analogCh),1,jj);
    hold on
    plot(xt(1:1000:end), d(1:1000:end));
    plot(xt([1 end]), [thr thr],'r','LineWidth',2);
    xlim([0 xt(end)]);
    ax = axis;
    if ~isempty(locsA)
        plot(locsA, ax(4),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none','MarkerSize',3);
    end
    if ~isempty(eventGroup)
        for kk = 1:size(eventGroup,1)
            plot([eventGroup(kk,1) eventGroup(kk,2)],[thr+100 thr+100],'LineWidth',10);
        end
    end
    xlabel('s'); ylabel(['Ch', num2str(analogCh(jj)),' (au)']); 
end
mkdir('Pulses');
saveas(gca,['pulses\analogPulsesDetection.png']);
% close all

filetarget = split(pwd,filesep); filetarget = filetarget{end};
if ~isempty(pul) % if no pulses, not save anything... 
    pulses.timestamps = stackCell(pul);
    pulses.amplitude = stackCell(val);
    pulses.duration = stackCell(dur);
    pulses.eventGroupID = stackCell(eventGroupID);
    pulses.analogChannel = stackCell(eventChannel);
    intsPeriods = [];
    for ii = 1:length(stimPer)
        intsPeriods = [intsPeriods; stimPer{ii}];
    end
    pulses.intsPeriods = intsPeriods;
    
    % sorting output
    [~, idx] = sort(pulses.timestamps(:,1));
    pulses.timestamps = pulses.timestamps(idx,:);
    pulses.amplitude = pulses.amplitude(idx,:);
    pulses.duration = pulses.duration(idx,:);
    pulses.eventGroupID = pulses.eventGroupID(idx,:);
    pulses.analogChannel = pulses.analogChannel(idx,:);
    
    [~, idx] = sort(pulses.intsPeriods(:,1));
    pulses.intsPeriods = pulses.intsPeriods(idx,:);
    
    % remove pulses that are too short (noise)
    if minDur
        ind = find(pulses.duration > minDur);
        
        pulses.timestamps = pulses.timestamps(ind,:);
        pulses.amplitude = pulses.amplitude(ind,:);
        pulses.duration = pulses.duration(ind,:);
        pulses.eventGroupID = pulses.eventGroupID(ind,:);
        pulses.analogChannel = pulses.analogChannel(ind,:);
        % need to deal with intsPeriods
    end
    
    disp('Saving locally...');
    save([filetarget '.pulses.events.mat'],'pulses');
else
    pulses = [];
end
 
cd(prevPath);
end

function [outMat] = stackCell(inCell)
    outMat = [];
    for ii = 1:length(inCell)
        outMat = [outMat; inCell{ii}'];
    end
end
