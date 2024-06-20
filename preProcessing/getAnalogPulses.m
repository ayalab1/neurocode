function [pulses] = getAnalogPulses(varargin)
%
% [[pul, val, dur] = getAnalogPulses(varargin)]
%
% [Find square pulses. If not argument it provide, it tries to find pulses
% in intan analog-in file]
%
%  INPUT
%
%    [parser]      [input parser, see below]
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values

% [forceDetect] [true or false to force detection (avoid load previous
%               detection, default false)]
% [analogCh]    [List of analog channels with pulses to be detected (it
%               support Intan Buzsaki Edition)]
% [data]        [R x C matrix with analog data. C is data, R should be
%               greater than 1]
% [samplingRate][Sampling frequency (in Hz), default 20000]
% [offset]      [Offset subtracted (in seconds), default 0]
% [periodLag]   [How long a pulse has to be far from other pulses to be
%               consider a different stimulation period (in seconds, default 20s)]
% [filename]    [File to get pulses from. Default, data file with folder
%               name in current directory]
% [manualThr]   [Check manually threslhold amplitude (default, false)]
% [groupPulses] [Group manually train of pulses (default, false)]
% [basepath]    [Path with analog data files to get pulses from]
% [minDur]      [pulses with shorter duration than this are removed]
% [showFig]     [Whether or not to show final pulse figure, best to set to
%               false when manually scoring long sessions. Default, true]
% [xaxis_adj_ints] [change xaxis for IDing groupPulses[xmin xmax]. Default 
%                  is [0 xmax]]
% [sessEpochs]  [session epochs during which to detect pulses. Should be 
%               integer indices for MergePoints. Default [], all epochs].
%
%
%  OUTPUT
%                 [pulses - events struct with the following fields
% [timestamps]    C x 2  matrix with pulse times in seconds. First column of C
%                 are the beggining of the pulses, second column of C are the end of
%                 the pulses]
% [amplitude]     [values of the pulses with respect baleline (normalized as 0)]
% [duration]      [Duration of the pulses. Note that default samplingRate is 20000]
% [eventID]       [Numeric ID for classifying various event types (C X 1)]
% [eventIDlabels] [label for classifying various event types defined in eventID (cell array C X 1)]
% [intsPeriods]   [Stimulation periods, as defined by perioLag]
%
% [concatenate_
% SEE ALSO
%
% [Manu-BuzsakiLab] [2018]
% [Antonio FR] [2021-2022]
% [HLarsson] [2023-2024]
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

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
addParameter(p,'showFig',true,@islogical);
addParameter(p,'forceDetect',false,@islogical);
addParameter(p,'xaxis_adj_ints', [], @isnumeric);
addParameter(p,'sessEpochs',[],@isnumeric);
addParameter(p,'smooth',0,@isnumeric);
addParameter(p,'maxValue',[],@isnumeric);



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
showFig = p.Results.showFig;
forceDetect = p.Results.forceDetect;
xaxis_adj_ints = p.Results.xaxis_adj_ints;
sessEpochs = p.Results.sessEpochs;
smooth = p.Results.smooth;
maxValue = p.Results.maxValue;


prevPath = pwd;
cd(basepath);

%%
filetarget = split(basepath,filesep); filetarget = filetarget{end};
if exist([filetarget '.pulses.events.mat'],'file') && ~forceDetect
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
        analogFile = split(basepath,'\'); analogFile = analogFile{end};
        analogFile = strcat(analogFile,'.dat');
    else
        analogFile = f.name;
    end
    
    if isempty(analogCh)
        error('No posible to run getAnalogPulses from Intan Buzsaki Ed with no analogCh inputs!');
    else
        analogCh = analogCh+1; % 0 to 1 index
    end
    parameters = LoadParameters(basepath); % read xml
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
        cd(localPaths(1).folder);
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

if manualThr && showFig
    h=figure;
end
% set(gcf,'Position',[100 -100 2500 1200]);
for jj = 1 : length(analogCh)
    fprintf(' ** Channel %3.i of %3.i... \n',jj, length(analogCh));
    disp('    Loading file...');
    if IntanBuzEd
        d = bz_LoadBinary(analogFile, 'frequency', samplingRate, 'nChannels', nChannels,'channels', analogCh(jj));
    else
        d = dataAnalogIn(analogCh(jj)-(min(analogCh)-1),:);
    end
    keyboard
    if ~isempty(maxValue), d(d>maxValue) = maxValue; end
    if smooth~=0, d = Smooth(double(d),smooth); end
    

    xt = linspace(0,length(d)/samplingRate,length(d));
    
    try % correct for different baselines in different subsessions (different rooms)
        MergePoints = getStruct(basepath,'MergePoints');
        for subsession = 1:size(MergePoints.timestamps,1)
            start = MergePoints.timestamps(subsession,1); stop = MergePoints.timestamps(subsession,2);
            in = xt>start & xt<stop;
            m = median(d(in)); % remove the median (baseline) signal
            d(in) = d(in) - m;
        end
    catch
        warning('error in attempting to remove subsession baseline: tell Raly');
    end
    
    % estimate if there are ANY pulses in the signal:
    % if there are pulses, we expect good separation between high signal (during pulse) and low signal (outside of pulse)
    % if there are no pulses, we expect more uniform / gaussian noise signal with poor separation (not bimodal)
    emThreshold = 0.3; % effectiveness metric goes from 0 (poor separation) to 1 (best separation); the number is an arbitrary threshold that Raly made up (having seen largest noise em-s of ~0.2, and lowest em-s of real pulses around >0.5)
    try % define "effectiveness metric" for the off-vs-on signal separation (pulse) for each subsession:
        MergePoints = getStruct(basepath,'MergePoints');
        for subsession = 1:size(MergePoints.timestamps,1)
            start = MergePoints.timestamps(subsession,1); stop = MergePoints.timestamps(subsession,2);
            in = xt>start & xt<stop;
            dd = Smooth(double(d(in)),smooth); midpoint = mean([min(dd) max(dd)]);
            em(subsession,1) = 1-(mean(dd>=midpoint)*var(dd(dd>=midpoint),1)+ mean(dd<midpoint)*var(dd(dd<midpoint),1)) / var(dd);
        end
    catch % no subsessions detected, use whole session:
        dd = double(d); midpoint = mean([min(dd) max(dd)]);
        em = 1-(mean(dd>=midpoint)*var(dd(dd>=midpoint),1)+ mean(dd<midpoint)*var(dd(dd<midpoint),1)) / var(dd);
        warning('error in attempting to remove subsession baseline: tell Raly');
    end
    
    if max(em)<emThreshold, % if ANY subsession has good pulses, that's enough for this condition to pass
        warning('Signal didn''t pass Raly''s threshold for detecting pulses. Aborting... Talk to Raly if you think this is a mistake and pulses are actually in the signal.');
        d(:) = 0; % set signal to zero to avoid tedious steps to find pulses
    end
    
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
        if ~isempty(sessEpochs)
            for sess_it = 1:size(sessEpochs,1)
                start = MergePoints.timestamps(sessEpochs(sess_it),1); stop = MergePoints.timestamps(sessEpochs(sess_it),2);
                hold on
                xline(start,'--g','LineWidth',2); xline(stop,'--r','LineWidth',2);
            end
        end
        xlabel('s'); ylabel('amp');
        title('Select threshold with the mouse and press left click...');
        drawnow
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
        [a, aa] = size(xaxis_adj_ints); %#ok<ASGLU>
        if aa ~=2
            error('Incorrect inputs for xlim! Change to [xmin xmax]')
        else
            disp('Adjusting xaxis for groupPulses!')
        end
        if aa == 2
            xlim(xaxis_adj_ints)
        else
            xlim([0 xt(end)])
        end
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
    
    try
        [~,temp(2,:)] = FindClosest(locsB,temp2,'higher'); % Raly: this is much faster in FindClosest is in the path
    catch
        for ii = 1 : length(temp2) % pair begining and end of the pulse
            try temp(2,ii) =  locsB(find(locsB - temp2(ii) ==...
                    min(locsB(locsB > temp2(ii)) - temp2(ii))));
            catch
                keyboard;
            end
        end
    end
    temp(:,find(temp(1,:) == 0)) = [];
    
    keyboard
    pul{jj} = gather(temp);
    clear temp temp2
    
    baseline_d = int32(median(d(1:100:end)));
    val{jj}=[];
    for ii = 1 : size(pul{jj},2) % value of the pulse respect basaline
        val{jj}(ii) = median(int32(d(int32(pul{jj}(1,ii) * samplingRate : pul{jj}(2,ii) * samplingRate)))) - baseline_d;
    end
    
    pul{jj} = pul{jj} - offset;
    % discard pulses < 2 * median(abs(x)/0.6745) as noise or pulses in negatives times
    idx = find((val{jj} < (thr-baseline_d)*0.4) | pul{jj}(1,:)<0);
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
    if manualThr && showFig
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
        if ~isempty(sessEpochs)
            for sess_it = 1:size(sessEpochs,1)
                start = MergePoints.timestamps(sessEpochs(sess_it),1); stop = MergePoints.timestamps(sessEpochs(sess_it),2);
                hold on
                xline(start,'--g','LineWidth',2); xline(stop,'--r','LineWidth',2);
            end
        end
    end
end
mkdir('Pulses');
saveas(gca,['pulses\analogPulsesDetection.png']);
% close all


filetarget = split(basepath,filesep); filetarget = filetarget{end};
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
    
    % remove pulses that are outside selected epochs
    if ~isempty(sessEpochs)
        ind = [];
        for i = 1:size(sessEpochs,1)
            ind1 = find((pulses.timestamps(:,1) >= MergePoints.timestamps(sessEpochs(i),1))&pulses.timestamps(:,2) <= MergePoints.timestamps(sessEpochs(i),2));
            ind = unique([ind; ind1]);
        end
        ind = sort(ind);
        pulses.timestamps = pulses.timestamps(ind,:);
        pulses.amplitude = pulses.amplitude(ind,:);
        pulses.duration = pulses.duration(ind,:);
        pulses.eventGroupID = pulses.eventGroupID(ind,:);
        pulses.analogChannel = pulses.analogChannel(ind,:);
    end
    
    disp('Saving locally...');
    save([filetarget '.pulses.events.mat'],'pulses');
else
    pulses = [];
end

cd(prevPath);


function [outMat] = stackCell(inCell)
outMat = [];
for ii = 1:length(inCell)
    outMat = [outMat; inCell{ii}'];
end




