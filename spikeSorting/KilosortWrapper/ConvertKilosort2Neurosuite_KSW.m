function ConvertKilosort2Neurosuite_KSW(rez)

% Converts KiloSort templates Klusta into klusters-compatible
% fet,res,clu,spk files.  Works on a single shank of a recording, assumes a
% 16bit .dat and an .xml file is present in "basepath" (home folder) and 
% that they are named basename.dat and basename.xml.  
% 
% Inputs:
%   basepath -  directory path to the main recording folder with .dat and .xml
%               as well as shank folders made by makeProbeMapKlusta2.m (default is
%               current directory matlab is pointed to)
%   basename -  shared file name of .dat and .xml (default is last part of
%               current directory path, ie most immediate folder name)

% Brendon Watson 2016
% Edited by Peter Petersen 2017

savepath = rez.ops.savepath;
basepath = rez.ops.basepath;
basename = rez.ops.basename;
if ~exist('basepath','var')
    [~,basename] = fileparts(cd);
    basepath = cd;
end
if ~exist('rez','var')
    load(fullfile(basepath,'rez.mat'))
end

Nchan = rez.ops.Nchan;
connected    = rez.connected;
xcoords      = rez.xc;
ycoords      = rez.yc;
% Nchan = rez.ops.Nchan;
% connected    = ones(Nchan, 1);
% xcoords      = ones(Nchan, 1);
% ycoords      = (1:Nchan)';

par = LoadXml(fullfile(basepath,[basename '.xml']));

totalch = par.nChannels;
sbefore = 16;%samples before/after for spike extraction
safter = 24;%... could read from SpkGroups in xml
if isfield(par.SpkGrps,'nSamples')
    if ~isempty(par.SpkGrps(1).nSamples);
        if isfield(par.SpkGrps,'PeakSample')
            if ~isempty(par.SpkGrps(1).PeakSample);
                sbefore = par.SpkGrps(1).PeakSample;
                safter = par.SpkGrps(1).nSamples - par.SpkGrps(1).PeakSample;
            end
        end
    end
end

if exist(rez.ops.fbinary,'file')
    datpath = rez.ops.fbinary;
else
    datpath = fullfile(basepath,[basename '.dat']);
end

% [spikeTimes, ii] = sort(spikeTimes);
spktimes = uint64(rez.st3(:,1));
clu = uint32(rez.st3(:,2));
amplitudes = rez.st3(:,3);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

mkdir(fullfile(savepath,'OriginalClus'))
%% do homework for assigning templates to shanks
% [~,shank]=fileparts(basepath);
templates = rez.Wraw;
% m = min(templates,[],2);%find the min value of each waveform on each channel
% [~,m] = min(m,[],1);%find which channel minimum is least
% m = squeeze(m);%which channel is minimum on each template
m = max(abs(templates),[],2);%find the most deviated value of each waveform on each channel
[~,m] = max(m,[],1);%find which channel has most deviated value for each templnate
m = squeeze(m);%squeeze to 1d vector

grouplookup = rez.ops.kcoords;%list of group/shank of each channel
templateshankassignments = grouplookup(m);%for the list of maximal channels, which group is each in 
allgroups = unique(grouplookup);

%Grp 0 contain discared channels
allgroups(allgroups==0) = [];

for groupidx = 1:length(allgroups)
    
    %if isfield(par.SpkGrps(groupidx),'Channels')
    %if ~isempty(par.SpkGrps(groupidx).Channels)
    % for each group loop through, find all templates clus
    tgroup          = allgroups(groupidx);%shank number
    ttemplateidxs   = find(templateshankassignments==tgroup);%which templates/clusters are in that shank
    ttemplates      = templates(:,:,ttemplateidxs);
    tPCFeatureInds  = pcFeatureInds(:,ttemplateidxs);
    
    tidx            = ismember(clu,ttemplateidxs);%find spikes indices in this shank
    tclu            = clu(tidx);%extract template/cluster assignments of spikes on this shank
    tspktimes       = spktimes(tidx);
    
    gidx            = find(rez.ops.kcoords == tgroup);%find all channels in this group
    channellist     = [];
    
    for ch = 1:length(par.SpkGrps)
        if ismember(gidx(1),par.SpkGrps(ch).Channels+1)
            channellist = par.SpkGrps(ch).Channels+1;
            break
        end
    end
    if isempty(channellist)
        disp(['Cannot find spkgroup for group ' num2str(groupidx) ])
        continue
    end
    
    %% spike extraction from dat
    if groupidx == 1;
        dat             = memmapfile(datpath,'Format','int16');
    end
    tsampsperwave   = (sbefore+safter);
    ngroupchans     = length(channellist);
    valsperwave     = tsampsperwave * ngroupchans;
    wvforms_all     = zeros(length(tspktimes)*tsampsperwave*ngroupchans,1,'int16');
    wvranges        = zeros(length(tspktimes),ngroupchans);
    wvpowers        = zeros(1,length(tspktimes));
    
    for j=1:length(tspktimes)
        try
            w = dat.data((double(tspktimes(j))-sbefore).*totalch+1:(double(tspktimes(j))+safter).*totalch);
            wvforms=reshape(w,totalch,[]);
            %select needed channels
            wvforms = wvforms(channellist,:);
    %         % detrend
    %         wvforms = floor(detrend(double(wvforms)));
            % median subtract
            wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
            wvforms = wvforms(:);
            
        catch
            disp(['Error extracting spike at sample ' int2str(double(tspktimes(j))) '. Saving as zeros']);
            disp(['Time range of that spike was: ' num2str(double(tspktimes(j))-sbefore) ' to ' num2str(double(tspktimes(j))+safter) ' samples'])
            wvforms = zeros(valsperwave,1);
        end

        %some processing for fet file
        wvaswv = reshape(wvforms,tsampsperwave,ngroupchans);
        wvranges(j,:) = range(wvaswv);
        wvpowers(j) = sum(sum(wvaswv.^2));

        lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
    %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,100000) == 0
            disp([num2str(j) ' out of ' num2str(length(tspktimes)) ' done'])
        end
    end
    wvranges = wvranges';
    
    %% Spike features
%     for each template, rearrange the channels to reflect the shank order
    tdx = [];
    for tn = 1:size(ttemplates,3);
        tTempPCOrder = tPCFeatureInds(:,tn);%channel sequence used for pc storage for this template
        for k = 1:length(channellist);
            i = find(tTempPCOrder==channellist(k));
            if ~isempty(i)
                tdx(tn,k) = i;
            else
                tdx(tn,k) = nan;
            end
        end
    end
    
    
    featuresperspike = 3; % kilosort default
    
    % initialize fet file
    fets    = zeros(sum(tidx),size(pcFeatures,2),ngroupchans);
    pct     = pcFeatures(tidx,:,:);
    
    %for each cluster/template id, grab at once all spikes in that group
    %and rearrange their features to match the shank order
    allshankclu = unique(tclu);
    
    for tc = 1:length(allshankclu)
        tsc     = allshankclu(tc);
        i       = find(tclu==tsc);
        tforig  = pct(i,:,:);%the subset of spikes with this clu ide
        tfnew   = tforig; %will overwrite
        
        ii      = tdx(tc,:);%handling nan cases where the template channel used was not in the shank
        gixs    = ~isnan(ii);%good vs bad channels... those shank channels that were vs were not found in template pc channels
        bixs    = isnan(ii);
        g       = ii(gixs);
        
        tfnew(:,:,gixs) = tforig(:,:,g);%replace ok elements
        tfnew(:,:,bixs) = 0;%zero out channels that are not on this shank
        try
            fets(i,:,:) = tfnew(:,:,1:length(par.SpkGrps(groupidx).Channels));
        catch
            keyboard
        end
    end
    %extract for relevant spikes only...
    % and heurstically on d3 only take fets for one channel for each original channel in shank... even though kilosort pulls 12 channels of fet data regardless
    tfet1 = squeeze(fets(:,1,1:length(par.SpkGrps(groupidx).Channels)));%lazy reshaping
    tfet2 = squeeze(fets(:,2,1:length(par.SpkGrps(groupidx).Channels)));
    tfet3 = squeeze(fets(:,3,1:length(par.SpkGrps(groupidx).Channels)));
    fets = cat(2,tfet1,tfet2,tfet3)';%     fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
%     fets = double(squeeze(fets(1,:,:)));
    %mean activity per spike
%     fetmeans = mean(fets,1);%this is pretty redundant with wvpowers
%     %find first pcs, make means of those... 
%     featuresperspike = 4;
%     firstpcslist = 1:featuresperspike:size(fets,1);
%     firstpcmeans = mean(fets(firstpcslist,:),1);
% 
%     nfets = size(fets,1)+1;
%     fets = cat(1,fets,fetmeans,firstpcmeans,wvpowers,wvranges,double(tspktimes'));
    fets = cat(1,double(fets),double(wvpowers),double(wvranges),double(tspktimes'));
    fets = fets';
    % fets = cat(1,nfets,fets);

    %% writing to clu, res, fet, spk

    cluname = fullfile(savepath, [basename '.clu.' num2str(tgroup)]);
    resname = fullfile(savepath, [basename '.res.' num2str(tgroup)]);
    fetname = fullfile(savepath, [basename '.fet.' num2str(tgroup)]);
    spkname = fullfile(savepath, [basename '.spk.' num2str(tgroup)]);
  %fet
    SaveFetIn(fetname,fets);

    %clu
    % if ~exist(cluname,'file')
        tclu = cat(1,length(unique(tclu)),double(tclu));
        fid=fopen(cluname,'w'); 
    %     fprintf(fid,'%d\n',clu);
        fprintf(fid,'%.0f\n',tclu);
        fclose(fid);
        clear fid
    % end

    %res
    fid=fopen(resname,'w'); 
    fprintf(fid,'%.0f\n',tspktimes);
    fclose(fid);
    clear fid

    %spk
    fid=fopen(spkname,'w'); 
    fwrite(fid,wvforms_all,'int16');
    fclose(fid);
    clear fid 

    disp(['Shank ' num2str(tgroup) ' done'])
    %end
    %end
end
clear dat
copyfile(fullfile(savepath, [basename,'.clu.*']),fullfile(savepath, 'OriginalClus'))

function SaveFetIn(FileName, Fet, BufSize);

if nargin<3 | isempty(BufSize)
    BufSize = inf;
end

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
  formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');
fprintf(outputfile, '%d\n', nFeatures);

if isinf(BufSize)
  
  temp = [round(100* Fet(:,1:end-1)) round(Fet(:,end))];
    fprintf(outputfile,formatstring,temp');
else
    nBuf = floor(size(Fet,1)/BufSize)+1;
    
    for i=1:nBuf 
        BufInd = [(i-1)*nBuf+1:min(i*nBuf,size(Fet,1))];
        temp = [round(100* Fet(BufInd,1:end-1)) round(Fet(BufInd,end))];
        fprintf(outputfile,formatstring,temp');
    end
end
fclose(outputfile);
