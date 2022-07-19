
function [optogeneticResponses] = getOptogeneticResponse(varargin)
% [optogeneticResponse] = getOptogeneticResponse(varargin)
%
% Computes Psth and a several statistical measures of the cell responses.
%
% <OPTIONALS>
% analogCh      List of analog channels with light pulses. If not provided,
%                   gets psth for all analog channels.
% digitalCh     List of digital channels with light pulses. By defaut,
%                   none.
% spikes        buzcode spikes structure, if not provided tries loadSpikes.
% basepath      By default pwd.
% numRep        For bootstrapping, default, 0 - no bootstrapping; recommended 500 if using.
% binSize       In seconds, default, 0.001.
% winSize       In seconds, default, 0.5.
% rasterPlot    Default true.
% ratePlot      Default true.
% winSizePlot   Default [-0.1 .5];
% plotStim      Default [0]. Can input stimulation onset and offset points
%               in seconds (ie [0 0.2]) to show stim times. If inputting an
%               empty array, nothing will plot.
%
% OUTPUTS
% optogeneticResponse
%
% Manu-BuzsakiLab 2021

% Parse options
p = inputParser;
addParameter(p,'analogCh',[],@isnumeric);
addParameter(p,'digitalCh',[],@isnumeric);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'numRep',0,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',1,@isnumeric);
addParameter(p,'rasterPlot',true,@islogical);
addParameter(p,'ratePlot',true,@islogical);
addParameter(p,'winSizePlot',[-.1 .5],@isnumeric);
addParameter(p,'plotStim',[0],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);

parse(p, varargin{:});
analogCh = p.Results.analogCh;
digitalCh = p.Results.digitalCh;
basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
rasterPlot = p.Results.rasterPlot;
ratePlot = p.Results.ratePlot;
saveMat = p.Results.saveMat;
winSizePlot = p.Results.winSizePlot;
plotStim = p.Results.plotStim;
force = p.Results.force;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.optogeneticResponse.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Optogenetic responses already computed! Loading file...');
    load(targetFile.name);
    return
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end
if isempty(analogCh)
    pulsesAnalog = getAnalogPulses;
else
    pulsesAnalog = getAnalogPulses('analogCh',analogCh);
end

pulsesDigital.timestamps = []; pulsesDigital.digitalChannel = [];
if ~isempty(digitalCh)
    session = getSession;
    digitalIn = getDigitalIn('all','fs',session.extracellular.sr); 
    for ii = 1:length(digitalCh)
        pulsesDigital.timestamps = [pulsesDigital.timestamps; digitalIn.ints{digitalCh(ii)}'];
        pulsesDigital.digitalChannel = [pulsesDigital.digitalChannel; ones(size(digitalIn.ints{digitalCh(ii)},2),1) * digitalCh(ii)];
    end
end

pulses.timestamps = [pulsesAnalog.timestamps; pulsesDigital.timestamps];  % combine pulses
pulses.channel = [pulsesAnalog.analogChannel; pulsesDigital.digitalChannel + max(pulsesAnalog.analogChannel)];  % combine pulses
pulses.analogChannel = [pulsesAnalog.analogChannel; nan(size(pulsesDigital.digitalChannel))];  % 
pulses.digitalChannel = [nan(size(pulsesAnalog.analogChannel)); pulsesDigital.digitalChannel];  % 
pulses.duration = round(pulses.timestamps(:,2) - pulses.timestamps(:,1),3);  % 

% get cell response
optogeneticResponses = [];
pulseDuration = median(pulses.duration);
discardNonFrequent = find(pulses.duration> pulseDuration+.001 | pulses.duration< pulseDuration-.001);
pulses.timestamps(discardNonFrequent,:) = [];
pulses.channel(discardNonFrequent,:) = [];

spikes = loadSpikes;
channels = unique(pulses.channel); % code per channel, channel x duration should be implemented... 
pulseDuration = min(pulseDuration); % because code only codes for channel, we take minimum duration channel for responses
timestamps_recording = min(pulses.timestamps(:,2)):1/1250:max(pulses.timestamps(:,2));
% generate random events for boostraping
disp('Generating boostrap template...');
nPulses = int32(length(pulses.channel));
randomEvents = [];
for kk = 1:numRep
    randomEvents{kk} = sort(randsample(timestamps_recording, nPulses))';
end
for ii = 1:length(spikes.UID)
    fprintf(' **Pulses from unit %3.i/ %3.i \n',ii, size(spikes.UID,2)); %\n
    if numRep > 0
        [stccg, t] = CCG([spikes.times{ii} randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate');
        t_duringPulse = t > 0 & t < pulseDuration; 
        randomRatesDuringPulse = nanmean(stccg(t_duringPulse,2:size(randomEvents,2)+1,1),1);
        optogeneticResponses.bootsTrapRate(ii) = mean(randomRatesDuringPulse);
        optogeneticResponses.bootsTrapRateStd(ii) = std(randomRatesDuringPulse);
        pd = fitdist(randomRatesDuringPulse','normal');
        optogeneticResponses.bootsTrapCI(ii,:) = pd.icdf([.001 0.999]);
    else
        optogeneticResponses.bootsTrapRate(ii) = NaN;
        optogeneticResponses.bootsTrapRateStd(ii) = NaN;
        optogeneticResponses.bootsTrapCI(ii,:) = [NaN NaN];
    end
    for jj = 1:length(channels)
        pul = pulses.timestamps(pulses.channel == channels(jj),:);
        if length(pul) > 100
            [stccg, t] = CCG({spikes.times{ii}, pul(:,1)},[],'binSize',binSize,'duration',winSize,'norm','rate');
            optogeneticResponses.responsecurve(ii,jj,:) = stccg(:,2,1);
            optogeneticResponses.responsecurveSmooth(ii,jj,:) = smooth(stccg(:,2,1));
            t_duringPulse = t > 0 & t < pulseDuration; 
            t_beforePulse = t > -pulseDuration & t < 0; 
            optogeneticResponses.responsecurveZ(ii,jj,:) = (stccg(:,2,1) - mean(stccg(t < 0,2,1)))/std(stccg(t < 0,2,1));
            optogeneticResponses.responsecurveZSmooth(ii,jj,:) = smooth((stccg(:,2,1) - mean(stccg(t < 0,2,1)))/std(stccg(t < 0,2,1)));
            optogeneticResponses.rateDuringPulse(ii,jj,1) = mean(stccg(t_duringPulse,2,1));
            optogeneticResponses.rateBeforePulse(ii,jj,1) = mean(stccg(t_beforePulse,2,1));
            optogeneticResponses.rateZDuringPulse(ii,jj,1) = mean(squeeze(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)));
%            [h, optogeneticResponses.modulationSignificanceLevel(ii,jj,1)] = kstest2(stccg(t_duringPulse,2,1),stccg(t_beforePulse,2,1));
            ci = optogeneticResponses.bootsTrapCI(ii,:);
            
            % Boostrap test
            if optogeneticResponses.rateDuringPulse(ii,jj,1) > ci(2)
                test = 1;
            elseif optogeneticResponses.rateDuringPulse(ii,jj,1) < ci(1)
                test = -1;
            else
                test = 0;
            end
            optogeneticResponses.bootsTrapTest(ii,jj,1) = test;
            
            % z-score change test
            if mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
                test = 1;
            elseif mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
                test = -1;
            else
                test = 0;
            end
            optogeneticResponses.zscoreTest(ii,jj,1) = test;
            
            % 3 ways test. If not boostrap, it would be 2 ways.
%             if (optogeneticResponses.rateDuringPulse(ii,jj,1) > ci(2) || isnan(ci(2))) && optogeneticResponses.modulationSignificanceLevel(ii,jj,1)<0.01...
%                     && mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
%                 test = 1;
%             elseif (optogeneticResponses.rateDuringPulse(ii,jj,1) < ci(1) || isnan(ci(1))) && optogeneticResponses.modulationSignificanceLevel(ii,jj,1)<0.01 ...
%                     && mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
%                 test = -1;
%             else
%                 test = 0;
%             end
             optogeneticResponses.threeWaysTest(ii,jj,1) = test;
             optogeneticResponses.channel(ii,jj,1) = channels(jj);
        else
            optogeneticResponses.responsecurve(ii,jj,:) = nan(duration/binSize + 1,1);
            optogeneticResponses.responsecurveZ(ii,jj,:) = nan(duration/binSize + 1,1);
            optogeneticResponses.modulationSignificanceLevel(ii,jj,1) = NaN;
            optogeneticResponses.rateDuringPulse(ii,jj,1) = NaN;
            optogeneticResponses.rateBeforePulse(ii,jj,1) = NaN;
            optogeneticResponses.rateZDuringPulse(ii,jj,1) = NaN;
            optogeneticResponses.bootsTrapTest(ii,jj,1) = NaN;
            optogeneticResponses.channel(ii,jj,1) = channels(jj);
        end
    end
    optogeneticResponses.timestamps = t;
end

optogeneticResponses.channel = mean(optogeneticResponses.channel,1);
optogeneticResponses.pulses = pulses;
optogeneticResponses.numRep = numRep;
optogeneticResponses.binSize = binSize;
optogeneticResponses.winSize = winSize;
optogeneticResponses.winSizePlot = winSizePlot;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.optogeneticResponse.cellinfo.mat'],'optogeneticResponses');
end

% PLOTS
% 1. Rasters plot
if rasterPlot
    t = optogeneticResponses.timestamps;
    for ii = 1:length(channels)
        st = pulses.timestamps(pulses.channel==channels(ii),1);
        if length(st) > 5000 % if more than 5000
            st = randsample(st, 5000);
            st = sort(st);
        end
        disp('   Plotting spikes raster and psth...');
        % [stccg, t] = CCG([spikes.times st],[],'binSize',0.005,'duration',1);
        figure;
        set(gcf,'Position',[200 -500 2500 1200]);
        for jj = 1:size(spikes.UID,2)
            fprintf(' **Pulses from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            rast_x = []; rast_y = [];
            for kk = 1:length(st)
                temp_rast = spikes.times{jj} - st(kk);
                temp_rast = temp_rast(temp_rast>winSizePlot(1) & temp_rast<winSizePlot(2));
                rast_x = [rast_x temp_rast'];
                rast_y = [rast_y kk*ones(size(temp_rast))'];
            end

            % spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
            resp = squeeze(optogeneticResponses.responsecurveSmooth(jj,ii,:));
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            plot(rast_x, rast_y,'.','MarkerSize',1)
            hold on
            plot(t(t>winSizePlot(1) & t<winSizePlot(2)), resp(t>winSizePlot(1) & t<winSizePlot(2)) * kk/max(resp)/2,'k','LineWidth',2); hold on
            if ~isempty(plotStim)
                if length(plotStim)==1
                    xline(plotStim,'r--');
                elseif length(plotStim)>2
                    error('length of plotStim exceeds length of 2. Please only define relative stimulation onset and offset times');
                else
                    xline(plotStim(1),'r--'); hold on
                    xline(plotStim(2),'r--');
                end
            end                    
            xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 kk]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('Trial');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (s)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\OptogenticRespRaster_ch',num2str(channels(ii)) ,'ch.png']); 
    end
end
% 2. Rate plot
if ratePlot
    t = optogeneticResponses.timestamps;
    for ii = 1:length(channels)
        figure
        subplot(1,2,1)
        imagesc([t(1) t(end)],[1 size(optogeneticResponses.responsecurve,1)],...
            squeeze(optogeneticResponses.responsecurveSmooth(:,optogeneticResponses.channel==channels(ii),:))); hold on 
        caxis([0 10]); colormap(jet);
        xline(0,'w--');
        set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
        title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
        
        subplot(1,2,2)
        imagesc([t(1) t(end)],[1 size(optogeneticResponses.responsecurve,1)],...
            squeeze(optogeneticResponses.responsecurveZSmooth(:,optogeneticResponses.channel==channels(ii),:))); hold on
        caxis([-3 3]); colormap(jet);
        if ~isempty(plotStim)
            if length(plotStim)==1
                xline(plotStim,'k--');
            elseif length(plotStim)>2
                error('length of plotStim exceeds length of 2. Please only define relative stimulation onset and offset times');
            else
                xline(plotStim(1),'k--'); hold on
                xline(plotStim(2),'k--');
            end
        end
        set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
        title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
        saveas(gcf,['SummaryFigures\AnalogPulsesPsth_',num2str(channels(ii)) ,'ch.png']); 
    end
end

cd(prevPath);
end