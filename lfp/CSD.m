
function [ csd ] = CSD (lfp, varargin)

% [ CSD ] = CSD (lfp, varargin)
% Calculates the 1D approximation of current source density (CSD) from a
% linear array of LFPs

% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       channels    vector with channels to inlcude. If empty take all (default)
%       win         time interval to compute CSD. If empty take all (default)
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       plotCSD     true/false. Default true.
%       plotLFP     true/false. Default true.
%    =========================================================================

% OUTPUT:
%    CSD           a buzcode structure with fields csd.data,
%                                                   csd.timestamps
%                                                   csd.samplingRate
%                                                   csd.channels 
%                                                   csd.params

% Antonio FR, 7/18

%% Parse inputs

p = inputParser;
addParameter(p,'channels',1:size(lfp.data,2),@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'win',[lfp.timestamps(1) lfp.timestamps(end)],@isnumeric);
addParameter(p,'spat_sm',0,@isnumeric);
addParameter(p,'temp_sm',0,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'rearrange_by_xml',true,@islogical);

parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;
basepath = p.Results.basepath;
rearrange_by_xml = p.Results.rearrange_by_xml;
%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

win = (p.Results.win*samplingRate)+1;


%% Compute CSD

lfp_frag = data(win(1):win(2),channels)*-1;

if rearrange_by_xml
    basename = basenameFromBasepath(basepath);
    % load channel map from basename.session
    session = loadSession(basepath,basename);
    lfp_frag = lfp_frag(:,[session.extracellular.electrodeGroups.channels{:}]); 
    % remove bad channels 
    lfp_frag(:,session.channelTags.Bad.channels) = []; 
    % save channels in mapped order and remove those excluded 
    channels = [session.extracellular.electrodeGroups.channels{:}];
    channels(:,session.channelTags.Bad.channels) = [];
end

% detrend
if doDetrend
   lfp_frag = detrend(lfp_frag')';
end
    
% temporal smoothing
if temp_sm > 0
   for ch = 1:size(lfp_frag,2) 
       lfp_frag(:,ch) = smooth(double(lfp_frag(:,ch)),temp_sm,'sgolay');
   end
end

% spatial smoothing
if spat_sm > 0
   for t = 1:size(lfp_frag,1) 
       lfp_frag(t,:) = smooth(double(lfp_frag(t,:)),spat_sm,'lowess');
   end
end

% calculate CSD 
CSD = diff(lfp_frag,2,2);

% generate output structure
csd.data = CSD;
csd.timestamps = timestamps(win(1):win(2));
csd.samplingRate = samplingRate;
csd.channels = channels; 
csd.params.spat_sm = spat_sm;
csd.params.temp_sm = temp_sm;
csd.params.detrend = doDetrend;

%% Plot

if plotLFP
    
    cmax = max(max(CSD)); 

    figure;
    subplot(1,2,1);
    contourf(timestamps(win(1):win(2)),1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD'); 
   
    subplot(1,2,2);
    for ch=1:size(lfp_frag,2)
        offset = 500*(ch-1);
        sh_tmp = (lfp_frag(:,ch)) + offset;
        plot(timestamps(win(1):win(2)),sh_tmp,'k','LineWidth',1.5); hold on;
        clear sh_tmp
    end
     set(gca,'YDir','reverse','YTickLabel',[]);
     axis tight
     xlim([csd.timestamps(1) csd.timestamps(end)]);
    xlabel('time (s)');ylabel('channel');title('LFP');   
    
elseif plotCSD  
    
     cmax = max(max(CSD)); 
   
     figure;
     contourf(timestamps(win(1):win(2)),1:size(CSD,2),CSD',40,'LineColor','none');hold on;
     colormap jet; caxis([-cmax cmax]);
     set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title(CSD); 
   
end

end


    




