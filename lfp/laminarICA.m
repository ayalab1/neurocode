function [ica] = laminarICA(varargin)
%    [ica] = doICA(varargin)
%
% Performs Independent Component Analysis (ICA) decomposition

%performs ICA of a laminar LFP profile using the logistic infomax ICA algorithm of Bell & Sejnowski (1995) with
% the natural gradient feature of Amari, Cichocki & Yang.
% The resulting indenpendent components (ICs)correspond to different
% synaptic sources of the LFPs (e.g. CA3 and entorhinal inputs for CA1
% LFPs). For more details about the method implementation see Fernandez-Ruiz,JNeurosci, 2012 and
% Schomburg et al., Neuron, 2014
%
% This function is still work in progress
%
% INPUTS
% <optional>
% basepath      Default pwd
% lfp           a  structure with fields   lfp.data,
%                                          lfp.timestamps
%                                          lfp.samplingRate
%                                          lfp.channels.
%               If not provided, runs getLFP('all') on basepath
%               IMPORTANT: lfp provided must be in the correct order of
%               channels! If not provided, this function will by default
%               load the channels in the correct order from the first
%               shank.
% shankNum      Specifies which shank to load channels from if lfp is not provided.
%               Default 1
% passband      Prefiltering passband interval, default [30 200]
% nICs          Number of independent components to extract. Default 8.
% saveMat       Save results, default true.
% force         Force analysis (disable loading option if already computed,
%                   default false)
% plotWeights   Default true
% plotCFC       Will also calculate theta-gamma CFC for each IC. Default true.
% regionChan    Optional, if a .hippocampalLayers.channelinfo.mat file
%               exists, or if input is provided, it will label these
%               channels on the ica output. Default order is or, pyr, rad,
%               slm
% chanRange     Optional. Specify whether you only want a subset of
%               channels from the shank to be used.
%
% TO DO: INCLUDE IMPORTANTS IMPUTS TO runica.m as additional arguments!
%
% OUTPUT
% ica           a buzcode structure with the following fields:
% .activations  independent components (this is the data matrix to use for
%               analysis)
% .data         independent components organized by explanatory variance.
% .timestamps
% .topo         Weight matrix
% .sphere       data sphering matrix (chans,chans)
% .weights      ICA weight matrix (comps,chans)
% .meanvar      Explained variance.
%
% Copyright (C) 2021 by Antonio FR. (using previous code from Manu V, Ipshita Z)
% This functions is a  wrapper for the EEGLAB function 'runica.mat' from Scott
% Makeig (CNL/The Salk Institute, La Jolla, 1996-). Copyright (C) 2004-2011
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Parse inputs
p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'shankNum',1,@isnumeric)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'passband',[30 200],@isnumeric)
addParameter(p,'nICs',15,@isnumeric)
addParameter(p,'plotWeights',true,@islogical)
addParameter(p,'plotCFC',true,@islogical)
addParameter(p,'regionChan',[],@isnumeric)
addParameter(p,'chanRange',[], @isnumeric)
addParameter(p,'thetaChannel',[],@isnumeric)

parse(p,varargin{:});
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
force = p.Results.force;
passband = p.Results.passband;
nICs = p.Results.nICs;
shankNum = p.Results.shankNum;
regionChan = p.Results.regionChan;
plotWeights = p.Results.plotWeights;
plotCFC = p.Results.plotCFC;
chanRange = p.Results.chanRange;
thetaChannel = p.Results.thetaChannel;

% Deal with inputs
prevBasepath = pwd;
cd(basepath); % will change this once I figure out which downstream functions require PWD

basename = basenameFromBasepath(basepath);

targetFile = dir([basepath,filesep,'*.ica.channelInfo.mat']);
if ~isempty(targetFile) && ~force
    disp('ICA already computed! Loading file.');
    load(targetFile.name);
    return
end

if isempty(lfp)
    try
        % load session to get anatomical groups (should be mapped)
        load([basepath,filesep,[basename,'.session.mat']])
        if isempty(chanRange)
            channelOrder =  session.extracellular.electrodeGroups(shankNum).channels{:};
        else
            channelOrder = session.extracellular.electrodeGroups(shankNum).channels{:}(chanRange);
        end
        
        % remove bad channels
        channelOrder(ismember(channelOrder,session.channelTags.Bad.channels(:))) = [];
        disp(['removing bad channels: ',num2str(session.channelTags.Bad.channels(:)')])
        % load lfp using mapping
        [lfp,infoLFP] = getLFP(channelOrder);
    catch
        error('LFP not found!');
    end
end

if isempty(regionChan)
    regFile = dir([basepath,filesep,'*.hippocampalLayers.channelinfo.mat']);
    if ~isempty(regFile)
        load(regFile.name);
        regionChan = hippocampalLayers.all;
    end
end

%% It is recommmend to filter LFPs in the band of interest (e.g. gamma)
disp('Filtering...');
lfpFilt = bz_Filter(lfp,'passband',passband,'order',4,'filter','butter'); %!!!

%% Run runica

[weights,sphere,meanvar,bias,signs,lrates,data] = runica(lfp(:,2:end)','lrate',1.0000e-03,'pca',nICs);

% Normalize sphere
sphere = sphere./norm(sphere);

% Generate mixing and unmixing matrix
unmixing = weights*sphere;
if (size(unmixing,1)==size(unmixing,2)) && rank(unmixing)==size(unmixing,1)
    mixing = inv(unmixing);
else
    mixing = pinv(unmixing);
end
activations = unmixing * lfp(:,2:end)';

%% Populate output structure
ica.activations = activations;
ica.data = data';
ica.timestamps = lfp(:,1);
ica.sphere = sphere;
ica.weights = weights;
ica.meanvar = meanvar;
ica.topo = mixing;
ica.unmixing = unmixing;
ica.samplingRate = infoLFP.samplingRate;
ica.channels = infoLFP.channels;

if saveMat
    disp('Saving results...');
    save([basepath,filesep,basename '.ica.channelInfo.mat'],'ica');
end

if plotWeights
    figure;
    plot(ica.topo,repmat([1:numel(ica.channels)]',1,nICs),'LineWidth',2);
    set(gca,'YDir','reverse')
    hold on
    legend({num2str(ica.channels)})
    ylim([1 size(ica.weights,2)])
    if ~isempty(regionChan)
        for ii = 1:length(regionChan)
            idx = find(ica.channels == regionChan(ii));
            if ~isempty(idx)
                line([-300 300],[idx idx],'Color',[0.5 0.5 0.5],'LineStyle','--');
            end
        end
    end
    if ~isfolder('ICA')
        mkdir('ICA')
    end
    saveas(gcf,[basepath,filesep,'ICA\Weights.png']);
end

if plotCFC
    if exist('hippocampalLayers','var')
        pyrCh = hippocampalLayers.pyramidal;
    else
        %Pick the middle channel from the shank for LFP
        pyrCh = channelOrder(floor(length(channelOrder)/2));
    end
    lfpTheta = getLFP(pyrCh,'interval',[500 1500]); % Only take 1000 seconds worth of data to keep the computation quick
    in = ica.timestamps>= 500 & ica.timestamps<= 1500;
    %% For each ICA, run CFC
    phaserange = 5:0.5:12;
    amprange = 30:5:200;
    % Make lfp structure, with first channel as the Pyr lfp, and remaining
    % channels as the ica components.
        
    lfpICA.timestamps = lfpTheta(:,1);
    lfpICA.data = lfpTheta(:,2:end);
    lfpICA.data(:,2:(size(ica.weights,1)+1)) = ica.activations(:,in)';
    lfpICA.channels = 1:1:(size(ica.weights,1)+1);
    [comodulogram] = CFCPhaseAmp(lfpICA,phaserange,amprange,'phaseCh',1,'ampCh',2:length(lfpICA.channels));
    
    set(gcf,'Position',[100 100 1400 300])
    saveas(gcf,[basepath,filesep,'ICA\CFC.png']);
end

cd(prevBasepath);
end
