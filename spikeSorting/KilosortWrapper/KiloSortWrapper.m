function savepath = KiloSortWrapper(varargin)
%
% Creates channel map from Neuroscope xml files, runs KiloSort and
% writes output data to Neurosuite format or Phy.
% 
% USAGE
%
% KiloSortWrapper
% Run from data folder. File basenames must be the
% same as the name as current folder
%
% KiloSortWrapper(varargin)
% Check varargin description below when input parameters are parsed
%
% Dependencies:  KiloSort (https://github.com/cortex-lab/KiloSort)
% 
% Copyright (C) 2016-2022 Brendon Watson and the Buzsakilab and Ralitsa Todorova (mahal threshold)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

disp('Running Kilosort spike sorting with the Buzsaki lab wrapper')

%% Parsing inputs
p = inputParser;
basepath = cd;
[~,basename] = fileparts(basepath);

addParameter(p,'basepath',basepath,@ischar)         % path to the folder containing the data
addParameter(p,'basename',basename,@ischar)         % file basenames (of the dat and xml files)
addParameter(p,'GPU_id',1,@isnumeric)               % Specify the GPU_id
addParameter(p,'SSD_path','D:\KiloSort',@ischar)    % Path to SSD disk. Make it empty to disable SSD
addParameter(p,'CreateSubdirectory',1,@isnumeric)   % Puts the Kilosort output into a subfolder
addParameter(p,'performAutoCluster',0,@isnumeric)   % Performs PhyAutoCluster once Kilosort is complete when exporting to Phy
addParameter(p,'mahal',12,@isnumeric)               % Perform mahalanobis threshold by default (set to Inf to not impose threshold)
addParameter(p,'config','',@ischar)                 % Specify a configuration file to use from the ConfigurationFiles folder. e.g. 'Omid'

parse(p,varargin{:})

basepath = p.Results.basepath;
basename = p.Results.basename;
GPU_id = p.Results.GPU_id;
SSD_path = p.Results.SSD_path;
CreateSubdirectory = p.Results.CreateSubdirectory;
performAutoCluster = p.Results.performAutoCluster;
config = p.Results.config;
mahalThreshold = p.Results.mahal;

if mahalThreshold<Inf, disp(['Will remove spikes exceeding a Mahalanobis distance of ' num2str(mahalThreshold)]); end

cd(basepath)

%% Checking if dat and xml files exist
if ~exist(fullfile(basepath,[basename,'.xml']))
    warning('KilosortWrapper  %s.xml file not in path %s',basename,basepath);
    return
elseif ~exist(fullfile(basepath,[basename,'.dat']))
    warning('KilosortWrapper  %s.dat file not in path %s',basename,basepath)
    return
end

%% Creates a channel map file
disp('Creating ChannelMapFile')
createChannelMapFile_KSW(basepath,basename,'staggered');

%% Loading configurations
%K%XMLFilePath = fullfile(basepath, [basename '.xml']);

if isempty(config)
    disp('Running Kilosort with standard settings')
    ops = KilosortConfiguration('basepath',basepath);
else
    disp('Running Kilosort with user specific settings')
    config_string = str2func(['KiloSortConfiguration_' config_version]);
    ops = config_string(XMLFilePath);
    clear config_string;
end

%% % Checks SSD location for sufficient space
if isfolder(SSD_path)
    FileObj = java.io.File(SSD_path);
    free_bytes = FileObj.getFreeSpace;
    dat_file = dir(fullfile(basepath,[basename,'.dat']));
    if dat_file.bytes*1.1<FileObj.getFreeSpace
        disp('Creating a temporary dat file on the SSD drive')
        ops.fproc = fullfile(SSD_path, [basename,'_temp_wh.dat']);
    else
        warning('Not sufficient space on SSD drive. Creating local dat file instead')
        ops.fproc = fullfile(basepath,'temp_wh.dat');
    end
else
    ops.fproc = fullfile(basepath,'temp_wh.dat');
end

%%
if ops.GPU
    disp('Initializing GPU')
    gpudev = gpuDevice(GPU_id); % initialize GPU (will erase any existing GPU arrays)
end
if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end

%% Lauches KiloSort
disp('Running Kilosort pipeline')
disp('PreprocessingData')
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization

disp('Fitting templates')
rez = fitTemplates(rez, DATA, uproj);  % fit templates iteratively

disp('Extracting final spike times')
rez = fullMPMU(rez, DATA); % extract final spike times (overlapping extraction)

%% Removes spikes that are beyond the threshold Mahalanobis distance before exporting to phy
if mahalThreshold<Inf
    nClusters = rez.ops.Nfilt;
    pcs = double(reshape(rez.cProjPC,size(rez.cProjPC,1),[]));
    nPCs = size(pcs,2);
    for i=1:nClusters
        ok = find(rez.st3(:,2)==i);
        if length(ok)>nPCs
            d = sqrt(mahal(pcs(ok,:),pcs(ok,:)));
            bad = ok(d>mahalThreshold);
            rez.st3(bad,2) = 0; % set bad spikes to unused cluster "0"
        end
    end
    rez.raw.st3 = rez.st3; rez.raw.cProj = rez.cProj; rez.raw.cProjPC = rez.cProjPC;  % save raw data from before these spikes were deleted

    % Remove the bad cluster
    bad = rez.st3(:,2)==0;
    disp(['Removing a total of ' num2str(sum(bad)) ' bad spikes (Mahalanobis threshold)']);
    rez.st3(bad,:) = [];
    rez.cProj(bad,:) = [];
    rez.cProjPC(bad,:) = [];
end

%% posthoc merge templates (under construction)
% save matlab results file

if CreateSubdirectory
    timestamp = ['Kilosort_' datestr(clock,'yyyy-mm-dd_HHMMSS')];
    savepath = fullfile(basepath, timestamp);
    mkdir(savepath);
    copyfile([basename '.xml'],savepath);
else
    savepath = fullfile(basepath);
end
rez.ops.basepath = basepath;
rez.ops.basename = basename;
rez.ops.savepath = savepath;
disp('Saving rez file')
% rez = merge_posthoc2(rez);
save(fullfile(savepath,  'rez.mat'), 'rez', '-v7.3');

%% export python results file for Phy
if ops.export.phy
    disp('Converting to Phy format')
    rezToPhy_KSW(rez);
    
%     % AutoClustering the Phy output
%     if performAutoCluster
%         PhyAutoClustering_km(savepath);
%     end
end

%% export Neurosuite files
if ops.export.neurosuite
    disp('Converting to Klusters format')
    load('rez.mat')
    rez.ops.root = pwd;
    clustering_path = pwd;
    basename = rez.ops.basename;
    rez.ops.fbinary = fullfile(pwd, [basename,'.dat']);
    Kilosort2Neurosuite(rez)
    
    writeNPY(rez.ops.kcoords, fullfile(clustering_path, 'channel_shanks.npy'));

    phy_export_units(clustering_path,basename);
end

%% Remove temporary file and resetting GPU
delete(ops.fproc);
reset(gpudev)
gpuDevice([])
disp('Kilosort Processing complete')

