function [events] = removeArtifactsFromEvents(events,varargin)
%removeArtifactsFromEvents - Remove artifacts from event according to the SD
% 
% INPUTS
%    'events'           Buzcode format events (i.e. ripples) structure.
%
% <optional>
%    'basepath'         Default 'pwd'
%    'events'           Structure containing the statistical test results.
%    'winSize'          .5
%    'figOpt'           Default true
%    'method'           'Minima' or 'std'
%    'stdThreshold'
% 
% OUTPUS
%    'events'           Buzcode format events (i.e. ripples) structure
%                           after event/spiking thresholing 
%
% Manu Valero - BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'winSize',.5);
addParameter(p,'figOpt',true,@islogical);
addParameter(p,'method','std');
addParameter(p,'stdThreshold',10);


parse(p,varargin{:});
basepath = p.Results.basepath;
winSize = p.Results.winSize;
figOpt = p.Results.figOpt;
method = p.Results.method;
stdThreshold = p.Results.stdThreshold;

prevPath = pwd;
cd(basepath);

%
if isfield(events.detectorinfo,'detectionchannel')
    lfp = getLFP(events.detectorinfo.detectionchannel);
else
    lfp = getLFP(events.detectorinfo.detectionparms.Channels(1));
end

disp('Computing ripples std...');
stdEvents = [];
parfor ii = 1:length(events.peaks) % compute rms
    idx = lfp.timestamps > events.peaks(ii)-winSize & lfp.timestamps < events.peaks(ii)+winSize;
    stdEvents(ii) = std(single(lfp.data(idx,1)));
end

if strcmpi(method,'minima')
    [N, edges] = histcounts(stdEvents);
    centers = edges(2:end) - diff(edges)/2;
    [~, locs] = findpeaks(smooth([0 N]),'MinPeakDistance',5);
    cutpoint = centers(min(N(locs(1):locs(2)))==N);
    cutpoint = cutpoint(1);
else
    cutpoint = mean(stdEvents) + std(stdEvents) * stdThreshold;
end
   
% 
validEvents = find(stdEvents<=cutpoint);
fields = fieldnames(events);
matchLength = size(events.timestamps,1);
for i = 1:length(fields)
    currentField = fields{i};
    if size(events.(currentField),1)==matchLength
        events.(currentField)= events.(currentField)(validEvents,:);
    elseif size(events.(currentField),2)==matchLength
        events.(currentField) = evetns.(currentField)(:,validEvents);
    end
end
      
events.artifactsRemovalParameters.cutpoint = cutpoint;
events.artifactsRemovalParameters.winSize = winSize;
events.removed = ~ismember((1:length(stdEvents)),validEvents);

fprintf('Keeping %4.0f of %4.0f events \n',length(validEvents),length(stdEvents));

if figOpt
    figure
    hold on
    histogram(stdEvents);
    ax = axis;
    plot([cutpoint cutpoint],ax(3:4),'r');
    ylabel('Counts'); xlabel('Event Std [Intan amplitude]');
    title('Artifact detection threshold','FontWeight','normal','FontSize',11);
    try
    saveas(gcf,'SummaryFigures\ArtifactDetectionThreshold.png');
    catch 
    saveas(gcf,'ArtifactDetectionThreshold.png');  
    end
end

cd(prevPath);

end