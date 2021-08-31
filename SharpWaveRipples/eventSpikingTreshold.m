
function [events] = eventSpikingTreshold(events,varargin)
% Descriptive and mean/median difference analysis, with serveral plot
% options.
% 
% INPUTS
%    'events'           Buzcode format events (i.e. ripples) structure.
%
% <optional>
%    'basepath'         Default 'pwd'
%    'spikes'           Buzcode format spikes structure. If not provided runs loadSpikes.      
%    'events'           Structure containing the statistical test results.
%    'spikingThreshold' .5 
%    'winSize'          .5
%    'eventSize'        .01
%    'figOpt'           Default true
% 
% OUTPUS
%    'events'           Buzcode format events (i.e. ripples) structure
%                           after event/spiking thresholing 
%
% Manu Valero - BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'spikes',[]);
addParameter(p,'spikingThreshold',.5);
addParameter(p,'winSize',.5);
addParameter(p,'eventSize',.01);
addParameter(p,'figOpt',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
spikes = p.Results.spikes;
spikingThreshold = p.Results.spikingThreshold;
winSize = p.Results.winSize;
eventSize = p.Results.eventSize;
figOpt = p.Results.figOpt;

% 
if isempty(spikes)
    spikes = loadSpikes;
end
[spikemat] = bz_SpktToSpkmat(spikes, 'dt',0.01,'overlap',6);
sSpkMat = zscore(sum(spikemat.data,2)/size(spikemat.data,2));

toCheck = find(events.peaks+winSize<spikemat.timestamps(end));
sizeResponse = length(int32(1:winSize*2/(mean(diff(spikemat.timestamps)))-1));
eventPopResponse = zeros(length(toCheck),sizeResponse);
parfor ii = 1: length(toCheck)
    temp = sSpkMat(spikemat.timestamps>=events.peaks(ii)-winSize ...
        & spikemat.timestamps<=events.peaks(ii)+winSize);
    if length(temp) < sizeResponse
        eventPopResponse(ii,:) = zeros(1,sizeResponse);
    else
        eventPopResponse(ii,:) = temp(int32(1:sizeResponse));
    end
end

t_event = linspace(-winSize,winSize,size(eventPopResponse,2));
eventResponse = zeros(size(events.peaks));
eventResponse(toCheck) = mean(eventPopResponse(:,t_event>-eventSize & t_event<eventSize),2);
[~,idx] = sort(eventResponse(toCheck));

% 
validEvents = find(eventResponse>spikingThreshold);
n_rips = size(events.timestamps,1);
for f = fields(events)'
    if size(events.(f{1}),1) == n_rips
        events.(f{1}) = events.(f{1})(validEvents,:);
    end
end
events.eventSpikingParameters.spikingThreshold = spikingThreshold;
events.eventSpikingParameters.winSize = winSize;
events.eventSpikingParameters.eventSize = eventSize;

fprintf('Keeping %4.0f of %4.0f events \n',length(validEvents),length(eventResponse));


if figOpt
    figure
    subplot(1,3,[1 2])
    hold on
    imagesc(t_event,1:length(eventPopResponse),eventPopResponse(idx,:),[-3 3]); colormap(jet);
    plot([t_event([1 end])], length(find(eventResponse<spikingThreshold))* ones(2,1) ,'r','LineWidth',1); axis tight;
    ylabel('Events'); xlabel('Time (s)'); set(gca,'YDir','normal','TickDir','out');

    subplot(1,3,3)
    hold on
    plot(eventResponse(idx),1:length(eventPopResponse)); 
    plot([spikingThreshold spikingThreshold], [1 length(eventPopResponse)],'r');
    xlabel('Response (SD)'); ylim([1 length(eventPopResponse)]); set(gca,'YDir','normal','TickDir','out');
    try 
        saveas(gcf,'SummaryFigures\eventSpikingThreshold.png');
    end
end

end