function PSTH = computePSTH(event,spikes,varargin)
% This is a generalized way for creating a PSTH for units for various events
% 
% INPUTS
% event  : event times formatted according to the CellExplorer's convention
% spikes : spikes formatted according to the CellExplorer's convention
% 
% OUTPUT
% psth
% 
% Dependencies: CCG

% AntonioFR, 10/21. Based on Peter Petersen's calc_PSTH

p = inputParser;

addParameter(p,'binCount',100,@isnumeric);        % how many bins (for half the window)
addParameter(p,'alignment','onset',@ischar);    % alignment of time ['onset','center','peaks','offset']
addParameter(p,'binDistribution',[0.4,0.2,0.4],@isnumeric);  % How the bins should be distributed around the events, pre, during, post. Must sum to 1
addParameter(p,'duration',0.15,@isnumeric);        % duration of PSTH (for half the window - used in CCG) [in seconds]
addParameter(p,'smoothing',5,@isnumeric);       % any gaussian smoothing to apply? units of bins.
addParameter(p,'percentile',99,@isnumeric);     % if events does not have the same length, the event duration can be determined from percentile of the distribution of events
addParameter(p,'plots',true,@islogical);        % Show plots?
addParameter(p,'eventName','',@ischar);         % Title used for plots
addParameter(p,'maxWindow',10,@isnumeric);      % Maximum window size in seconds
addParameter(p,'zscorePlot',true,@islogical);   % plot z-scored response

parse(p,varargin{:})

binCount = p.Results.binCount;
alignment = p.Results.alignment;
binDistribution = p.Results.binDistribution;
duration = p.Results.duration;
smoothing = p.Results.smoothing;
percentile = p.Results.percentile;
eventName = p.Results.eventName;
plots = p.Results.plots;
maxWindow = p.Results.maxWindow;
zscorePlot = p.Results.zscorePlot;

% If no duration is given, an optimal duration is determined
if duration == 0
    durations = diff(event.timestamps');
    stim_duration = prctile(sort(durations),percentile);
    duration = min(max(round(stim_duration*1000),50)/1000,maxWindow);
end

binSize = max(round(duration/binCount*1000),1)/1000; % minimum binsize is 0.5ms.

% Determine event alignment
switch alignment
    case 'onset'
        event_times = event.timestamps(:,1);
        padding = binDistribution(1)/binDistribution(2)*duration;
        %binsToKeep = int64(ceil((padding+duration/2)/binSize):(duration+padding)*2/binSize);
         binsToKeep = int64(ceil(padding/binSize):ceil((duration*2+padding)/binSize));
    case 'center'
        event_times = mean(event.timestamps);
        padding = 0;
        binsToKeep = 1:duration*2/binSize;
    case 'offset'
        event_times = event.timestamps(:,2);
        padding = binDistribution(3)/binDistribution(2)*duration;
        %binsToKeep = 1:(duration+padding)*2/binSize-ceil((padding+duration/2)/binSize);
        binsToKeep = int64(ceil(padding/binSize):ceil((duration*2+padding)/binSize));
    case 'peaks'
        event_times = event.peaks;
        padding = 0;
        binsToKeep = 1:duration*2/binSize;
end

disp(['  ', num2str(length(event_times)), '  events, duration set to: ', num2str(duration), ' sec, aligned to ', alignment,', with binsize: ' num2str(binSize)])

% Determining the bins interval for metrics
binsPre = 1:floor(binDistribution(1)*length(binsToKeep));
binsEvents = floor(binDistribution(1)*length(binsToKeep))+1:floor((binDistribution(1)+binDistribution(2))*length(binsToKeep));
binsPost = floor((binDistribution(1)+binDistribution(2))*length(binsToKeep))+1:length(binsToKeep);

% Calculating PSTH
PSTH_out = [];
for j = 1:numel(spikes.times)
    [spike_times,index] = sort([spikes.times{j};event_times(:)]);
    spike_cluster_index = [ones(size(spikes.times{j}));2*ones(size(event_times(:)))];
    [ccg,time] = CCG(spike_times,spike_cluster_index(index),'binSize',binSize,'duration',(duration+padding)*2);
    PSTH_out(:,j) = ccg(binsToKeep+1,2,1)./numel(event_times)/binSize;
end
time = time(binsToKeep+1);

modulationIndex = mean(PSTH_out(binsEvents,:))./mean(PSTH_out(binsPre,:));
modulationSignificanceLevel = [];
for i = 1:size(PSTH_out,2)
    [~,p_kstest2] = kstest2(PSTH_out(binsEvents,i),PSTH_out(binsPre,i));
    modulationSignificanceLevel(i) = p_kstest2;
end

if smoothing>0
    PSTH_out = nanconv(PSTH_out,ce_gausswin(smoothing)/sum(ce_gausswin(smoothing)),'edge');
end

[~,modulationPeakResponseTime] = max(PSTH_out);
modulationPeakResponseTime = time(modulationPeakResponseTime);

PSTH.responsecurve = PSTH_out;
PSTH.time = time;
PSTH.alignment = alignment;

PSTH.modulationIndex = modulationIndex;
PSTH.modulationPeakResponseTime = modulationPeakResponseTime';
PSTH.modulationSignificanceLevel = modulationSignificanceLevel;

if plots
    %[~,index2] = sort(modulationIndex,'descend');
    [~,index3] = sort(modulationPeakResponseTime);    
    figure,
    subplot(2,1,1);
    plot(time,mean(PSTH_out')','LineWidth',2); hold on; 
    plot(time,mean(PSTH_out')+std(PSTH_out'),'--b');hold on; 
    plot(time,mean(PSTH_out')-std(PSTH_out'),'--b');hold on;
    xline(0,'--k');hold on; ylabel('mod. index');
    title(eventName)
    subplot(2,1,2)
    if zscorePlot
    imagesc(time,[1:size(PSTH_out,2)],zscore(PSTH_out(:,index3))',[-3 3]), xlabel('time'), ylabel('units');hold on;
    else
    imagesc(time,[1:size(PSTH_out,2)],(PSTH_out(:,index3))'), xlabel('time'), ylabel('units');hold on;        
    end
    xline(0,'--k');hold on; 
end
