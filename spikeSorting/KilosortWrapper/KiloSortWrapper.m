function savepath = KiloSortWrapper(varargin)
%
%  USAGE
%
%   KiloSortWrapper
%   Run from data folder. File basenames must be the
%   same as the name as current folder
%
%   KiloSortWrapper(varargin)
%   Check varargin description below when input parameters are parsed
%
%  INPUT
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'basepath'            basepath for wrapper (if not listed, cd)
%     'p'                   parser
%     '.basepath'           path to the folder containing the data
%     '.basename'           file basenames (of the dat and xml files)
%     '.GPU_id'             specify the GPU_id (default is 1)
%     '.rejectchannels'     specify list of channels to ignore while spike... 
%                           ...sorting (base 1, add 1 to neuroscope numbering)
%     '.SSD_path'           path to SSD disk. Make it empty to disable SSD
%     '.CreateSubdirectory' puts the Kilosort output into a subfolder...
%                           ...(defualt is 1)
%     '.performAutoCluster' performs PhyAutoCluster once Kilosort is... 
%                           ...complete when exporting to Phy 
%     '.config'             specify a configuration file to use from the...
%                           ....ConfigurationFiles folder. e.g. 'Omid'
%     '.NT'                 specify desired batch size...
%                           ...(default = 32*1024; reduce if out of memory)
%     =========================================================================
%
%  OUTPUT
%
%    'rez'                  rez structure for input into KiloSort
%  
%  NOTE
%
%  SEE
%
%
%   Dependencies:  KiloSort (https://github.com/cortex-lab/KiloSort)
%                  createChannelMapFile_KSW, KilosortConfiguration,
%                  convertOpenEphysToRawBInary, gpuDevice, preprocessData
%                  fitTemplates, fullMPMU, rezToPhy_KSW, Kilosort2Neurosuite
%                  writeNPY, phy_export_units
%
% % Copyright (C) 2016-2022 Brendon Watson and the Buzsakilab
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.Creates channel map from Neuroscope xml files, runs KiloSort and
% writes output data to Neurosuite format or Phy.
% 

disp('Running Kilosort spike sorting with the Buzsaki lab wrapper')

%% Parsing inputs
p = inputParser;
basepath = cd;
[~,basename] = fileparts(basepath);

addParameter(p,'basepath',basepath,@ischar)         % path to the folder containing the data
addParameter(p,'basename',basename,@ischar)         % file basenames (of the dat and xml files)
addParameter(p,'GPU_id',1,@isnumeric)               % Specify the GPU_id
addParameter(p,'rejectchannels',[],@isnumeric)      % Specify list of channels to ignore while spike sorting (base 1, add 1 to neuroscope numbering)
addParameter(p,'SSD_path','D:\KiloSort',@ischar)    % Path to SSD disk. Make it empty to disable SSD
addParameter(p,'CreateSubdirectory',1,@isnumeric)   % Puts the Kilosort output into a subfolder
addParameter(p,'performAutoCluster',0,@isnumeric)   % Performs PhyAutoCluster once Kilosort is complete when exporting to Phy
addParameter(p,'config','',@ischar)                 % Specify a configuration file to use from the ConfigurationFiles folder. e.g. 'Omid'
addParameter(p,'NT',[],@isnumeric)                  % Specify desired batch size (default = 32*1024; reduce if out of memory)

parse(p,varargin{:})

basepath = p.Results.basepath;
basename = p.Results.basename;
GPU_id = p.Results.GPU_id;
SSD_path = p.Results.SSD_path;
CreateSubdirectory = p.Results.CreateSubdirectory;
performAutoCluster = p.Results.performAutoCluster;
config = p.Results.config;
rejectChannels = p.Results.rejectchannels;
NT = p.Results.NT;

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
createChannelMapFile_KSW(basepath,basename,'staggered',rejectChannels);

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

if ~isempty(NT)
    ops.NT = NT + ops.ntbuff;
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

