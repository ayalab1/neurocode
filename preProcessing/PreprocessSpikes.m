function PreprocessSpikes(basepath, varargin)
%
%
% [PreprocessSpikes] - [secondary preprocessing for after manual
%                       clustering]
%
%  [Secondary preprocessing function for after KiloSort and manual
%  clustering has been completed. This processes the spike waveform and
%  runs initial cell metrics]
%
%  INPUTS
%    [basepath]     [Character input with the file path for your files
%                   (specifically session file)]
%    <options>     optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%   [multiKilosort]  [logical value denoting whether or not there are
%                    multiple KiloSort files to be concatenated together]
%   [showGUI]        [logical value for showing the session template]
%   [showCellMet]    [logical value for whether or not to pull up cell
%                    metrics GUI after running initial cell metrics]
%   [prePhy]         [logical value denoting whether or not this is a
%                    temporary preprocessing run. Default false]
%   [spikeLabels]    [cell structure of spike sorting labels to read in.
%                    Default is good. It is recommended to set to {'good',
%                    'unsorted'} when prePhy=true]
%    =========================================================================
%
%  OUTPUTS
%    N/A
%
% SEE ALSO
%
% [HeathLarsson] [2021-2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p,'multiKilosort',false,@islogical);
addParameter(p,'showGUI',false,@islogical);
addParameter(p,'showCellMet',true,@islogical);
addParameter(p,'prePhy',false,@islogical);
addParameter(p,'spikeLabels',{'good'},@iscell);
parse(p,varargin{:});

multiKilosort = p.Results.multiKilosort;
showGUI = p.Results.showGUI;
showCellMet = p.Results.showCellMet;
prePhy = p.Results.prePhy;
spikeLabels = p.Results.spikeLabels;

%% 1- extract spike times and waveforms for sorted clusters
cd(basepath);
basename = basenameFromBasepath(pwd);
%if you intend to concatenate multiple kilosort runs into one spikes structure
try 
    load(fullfile(basepath,[basename '.session.mat']),'session');
    session.extracellular.chanCoords; % if this already exists, no need to re-load rez file
catch
    session = sessionTemplate(basepath,'showGUI',showGUI);
    save(fullfile(basepath,[basename '.session.mat']),'session');
end
f = dir('Kilosort*');
if (size(f,1) ~= 1)&&(~multiKilosort)
    error('Too many kiloSort folders - should multiKilosort=1?');
elseif (size(f,1) ~= 1)&&(multiKilosort)
    warning('Make sure all spike groups in neuroscope are reactivated (yellow tab, groups should be assigned to a numbered group');
end
% Make sure there is only one KiloSort folder before running, unless you
% needed to spike sort probes separately (multiKilosort=1).

% The hippocampal KiloSort folder should be listed first in the session
% folder for organization, but this is not necessary for running

% check if spikes.cellinfo has already been created
pre_exist_spike_files = dir('*spikes*.cellinfo.mat');
if (size(f,1) > 1) && multiKilosort
    for i = 1:size(f,1)
        session.spikeSorting{1,1}.relativePath = f(i).name;
        % if all the spike files exist, load instead of making from scratch
        if length(pre_exist_spike_files) == length(f)
            temp_spikes = load(pre_exist_spike_files(i).name);
            tempSpikes{i} = temp_spikes.spikes;
        else
            tempSpikes{i} = loadSpikes('session',session,...
                'clusteringpath',[f(i).folder filesep f(i).name],...
                'savemat', false,'labelsToRead',spikeLabels);
        end
        if i == 1
            spikes = tempSpikes{i};
        else
            fields{1}= fieldnames(spikes);
            fields{2} = fieldnames(tempSpikes{i});
            if length(fields{1}) >= length(fields{2})
                base = 1; comp = 2;
            else
                base = 2; comp = 1;
            end
            prevUID = length(spikes.UID);
            for j = 1:length(fields{base})
                currentField = fields{base}{j};
                if sum(contains(fields{comp},currentField))
                    %check special cases
                    if contains(currentField, 'spindices')
                        tempSpindices2 = tempSpikes{i}.spindices(:,2)+prevUID;
                        spikes.spindices = [spikes.spindices; [tempSpikes{i}.spindices(:,1) tempSpindices2]];
                    elseif contains(currentField, 'UID')
                        spikes.UID = [spikes.UID tempSpikes{i}.UID+prevUID];
                    elseif contains(currentField, 'basename')
                        if ~(spikes.basename == tempSpikes{i}.basename)
                            error('Incompatible basenames');
                        end
                    elseif contains(currentField, 'numcells')
                        spikes.numcells = spikes.numcells+tempSpikes{i}.numcells;
                    elseif contains(currentField, 'sr')
                        if ~(spikes.sr == tempSpikes{i}.sr)
                            error('Incompatible sampling rates');
                        end
                    elseif contains(currentField, 'processinginfo')
                        disp('Processing info assumed to be the same')
                    elseif iscell(spikes.(currentField))
                        spikes.(currentField) = cat(2,spikes.(currentField),tempSpikes{i}.(currentField));
                    else
                        spikes.(currentField) = [spikes.(currentField) tempSpikes{i}.(currentField)];
                    end
                else
                    disp(['Forced to remove ' currentField]);
                    if isfield(spikes, currentField)
                        spikes = rmfield(spikes,currentField);
                    end
                end
            end
        end
    end
    if prePhy
        save(strcat(basename,'.unsorted.spikes.cellinfo.mat'),'spikes');
    else
        save(strcat(basename,'.spikes.cellinfo.mat'),'spikes');
    end
    clear tempSpikes i j tempSpindices2 prevUID base comp
elseif (size(f,1)==1)&&(multiKilosort)
    error('Only one kilosort folder present, cannot combine multiple runs');
else
    if prePhy
        spikes = loadSpikes('session',session,'clusteringpath',[f.name],'labelsToRead',spikeLabels,'basename',[basename '.unsorted'],'getWaveformsFromDat',false);
    else
        spikes = loadSpikes('session',session,'clusteringpath',[f.name],'labelsToRead',spikeLabels);
    end
end
%% 2 - compute basic cell metrics
if exist([basepath '\anatomical_map.csv'])
    channel_mapping;
end
if prePhy
    spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name],'labelsToRead',spikeLabels,'basename',[basename '.unsorted'],'getWaveformsFromDat',false,'forceReload',true);
    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial','monoSynaptic_connections'},'saveAs','unsorted.cell_metrics','getWaveformsFromDat',false);
else
    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,'manualAdjustMonoSyn',false);
end
% GUI to manually curate cell classification
if showCellMet
    cell_metrics = CellExplorer('metrics',cell_metrics);
end
end