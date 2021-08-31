function out = hackInfo(varargin)
% patch for dealing with different session info formats - just call this to
% get basic ephys info (could add more whenever)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr)
addParameter(p,'useOld',false,@islogical)


parse(p,varargin{:})
basepath = p.Results.basepath;
useOld = p.Results.useOld;


try
    checkFile('basepath',basepath,'fileType','.session.mat');
    newForm = true;
catch
    try
        chcekFile('basepath',basepath,'fileType','.sessionInfo.mat');
        newForm = false;
    catch
        error('Need to calculate some kind of session info file')
    end
end

if newForm && ~useOld
    session = getSession('basepath',basepath);
    out.OGformat = 'new';
    out.nChannel = session.extracellular.nChannels;
    out.one.channels = 1:out.nChannel;
    try
        out.one.badChannels = session.channelTags.Bad.channels;
    catch
        out.one.badChannels = [];
    end
    out.zero.channels = out.one.channels-1;
    out.zero.badChannels = out.one.badChannels-1;
    out.nShank = session.extracellular.nElectrodeGroups;
    out.lfpSR = session.extracellular.srLfp;
    out.spikeSR = session.extracellular.sr;
    out.one.AnatGrps = session.extracellular.electrodeGroups.channels;
    out.zero.AnatGrps = cellfun(@(x)x-1,out.one.AnatGrps,'UniformOutput',false);
else
    sessionInfo = bz_getSessionInfo(basepath);
    out.OGformat = 'old';
    out.nChannel = sessionInfo.nChannels;
    out.zero.channels = sessionInfo.channels;
    out.zero.badChannels = sessionInfo.badchannels;
    out.one.channels = out.zero.channels+1;
    out.one.badChannels = out.zero.badChannels+1;
    out.nShank = sessionInfo.nElecGps;
    out.lfpSR = sessionInfo.lfpSampleRate;
    out.spikeSR = sessionInfo.rates.wideband;
    for shIdx = 1:out.nShank
        out.zero.AnatGrps{shIdx} = sessionInfo.AnatGrps(shIdx).Channels;
    end
    out.one.AnatGrps = cellfun(@(x)x+1,out.zero.AnatGrps,'UniformOutput',false);
end

shankID = zeros(1,out.nChannel);
for shIdx = 1:out.nShank
    shankID(out.one.AnatGrps{shIdx}) = shIdx;
end

out.shankID = shankID;
    
    
end
    
    