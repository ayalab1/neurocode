
function [events] = removeArtifactsFromEvents(events,varargin)
% Remove artifacts from event according to the SD.
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
addParameter(p,'stdThreshold',2);


parse(p,varargin{:});
basepath = p.Results.basepath;
winSize = p.Results.winSize;
figOpt = p.Results.figOpt;
method = p.Results.method;
stdThreshold = p.Results.stdThreshold;

prevPath = pwd;
cd(basepath);

%
lfp = bz_GetLFP(events.detectorinfo.detectionchannel);

disp('Computing ripples std...');
stdEvents = [];
parfor ii = 1:length(events.peaks) % compute rms
%    rmsEvents(ii) = rms(double(lfp.data(find(lfp.timestamps>events.peaks(ii)-winSize &...
%        lfp.timestamps<events.peaks(ii)+winSize))));
    stdEvents(ii) = std(double(lfp.data(find(lfp.timestamps>events.peaks(ii)-winSize &...
        lfp.timestamps<events.peaks(ii)+winSize))));
end

if strcmpi(method,'minima')
    [N, edges] = histcounts(stdEvents);
    centers = edges(2:end) - diff(edges)/2;
    [~, locs] = findpeaks(smooth([0 N]),'MinPeakDistance',5);
    cutpoint = centers(find(min(N(locs(1):locs(2)))==N));
    cutpoint = cutpoint(1);
else
    cutpoint = mean(stdEvents) + std(stdEvents) * stdThreshold;
end
   
% 
validEvents = find(stdEvents<=cutpoint);
events.timestamps = events.timestamps(validEvents,:);
try 
    events.peaks = events.peaks(validEvents,:);
    events.peakNormedPower = events.peakNormedPower(validEvents,:);
end
events.artifactsRemovalParameters.cutpoint = cutpoint;
events.artifactsRemovalParameters.winSize = winSize;

fprintf('Keeping %4.0f of %4.0f events \n',length(validEvents),length(stdEvents));

if figOpt
    figure
    hold on
    histogram(stdEvents);
    ax = axis;
    plot([cutpoint cutpoint],ax(3:4),'r');
    ylabel('Counts'); xlabel('Event Std [Intan amplitude]');
    title('Artifact detection threshold','FontWeight','normal','FontSize',11);
    saveas(gcf,'SummaryFigures\ArtifactDetectionThreshold.png');
end

cd(prevPath);

end