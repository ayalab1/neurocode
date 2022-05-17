% function preprocessSpikes(session)
% This is a wrapper to concantenate the basic spike related functions,
% meant to be run right after the manual clustering

%% 0 - Spike sorting
% Automatic spike sorting should be run with Kilosort (using preprocessSession)
% After that, copy yo your local SSD the Kilosort folder, .dat and .xml files
% Do manual curation using Phy2

%% 1- extract spike times and waveforms for sorted clusters
ifMultiKiloSort = 0; %if you intend to concatenate multiple kilosort runs into one spikes structure
session = sessionTemplate(pwd,'showGUI',false);
f = dir('Kilosort*');
% Make sure there is only one KiloSort folder before running, unless you needed to spike sort probes separately
% The hippocampal KiloSort folder should be listed first in the session folder for organization, but this is not necessary for running
if (size(f,1) > 1) && ifMultiKiloSort
    for i = 2:size(f,1)
        session.spikeSorting{1,1}.relativePath = f(i).name;
        tempSpikes{i} = loadSpikes('session',session,'clusteringpath',[f(i).folder filesep f(i).name], 'savemat', false);
        if i == 1
            spikes = tempSpikes{i};
        else
            spikes.ids = cat(2,spikes.ids,tempSpikes{i}.ids);
            spikes.ts = cat(2,spikes.ts,tempSpikes{i}.ts);
            spikes.times = cat(2, spikes.times, tempSpikes{i}.times);
            spikes.cluID = [spikes.cluID tempSpikes{i}.cluID];
            spikes.total = [spikes.total tempSpikes{i}.total];
            spikes.amplitudes = cat(2,spikes.amplitudes,tempSpikes{i}.amplitudes);
            spikes.maxWaveformCh = [spikes.maxWaveformCh tempSpikes{i}.maxWaveformCh];
            spikes.maxWaveformCh1 = [spikes.maxWaveformCh1 tempSpikes{i}.maxWaveformCh1];
            if isfield(tempSpikes{i},'phy_maxWaveformCh1')
                spikes.phy_maxWaveformCh1 = [spikes.phy_maxWaveformCh1 tempSpikes{i}.phy_maxWaveformCh1];
            end
            spikes.phy_amp = [spikes.phy_amp tempSpikes{i}.phy_amp];
            spikes.numcells = spikes.numcells+tempSpikes{i}.numcells;
            tempSpindices2 = tempSpikes{i}.spindices(:,2)+length(spikes.UID);
            spikes.spindices = [spikes.spindices; [tempSpikes{i}.spindices(:,1) tempSpindices2]];
            spikes.UID = [spikes.UID tempSpikes{i}.UID+length(spikes.UID)];
            spikes.shankID = [spikes.shankID tempSpikes{i}.shankID];
            spikes.rawWaveform = cat(2,spikes.rawWaveform,tempSpikes{i}.rawWaveform);
            spikes.filtWaveform = cat(2,spikes.filtWaveform,tempSpikes{i}.filtWaveform);
            spikes.rawWaveform_all = cat(2,spikes.rawWaveform_all,tempSpikes{i}.rawWaveform_all);
            spikes.filtWaveform_all = cat(2,spikes.filtWaveform_all,tempSpikes{i}.filtWaveform_all);
            spikes.rawWaveform_std = cat(2,spikes.rawWaveform_std,tempSpikes{i}.rawWaveform_std);
            spikes.filtWaveform_std = cat(2,spikes.filtWaveform_std,tempSpikes{i}.filtWaveform_std);
            spikes.timeWaveform = cat(2,spikes.timeWaveform,tempSpikes{i}.timeWaveform);
            spikes.timeWaveform_all = cat(2,spikes.timeWaveform_all,tempSpikes{i}.timeWaveform_all);
            spikes.peakVoltage = [spikes.peakVoltage tempSpikes{i}.peakVoltage];
            spikes.channels_all = cat(2,spikes.channels_all, tempSpikes{i}.channels_all);
            spikes.peakVoltage_sorted = cat(2,spikes.peakVoltage_sorted,tempSpikes{i}.peakVoltage_sorted);
            spikes.maxWaveform_all = cat(2,spikes.maxWaveform_all,tempSpikes{i}.maxWaveform_all);
            spikes.peakVoltage_expFitLengthConstant = [spikes.peakVoltage_expFitLengthConstant tempSpikes{i}.peakVoltage_expFitLengthConstant];
        end
    end
    basename = basenameFromBasepath(pwd);
    save(strcat(basename,'.spikes.cellinfo.mat'),'spikes');
    clear tempSpikes i tempSpindices2
else
    spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
end
%% 2 - compute basic cell metrics
cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,'manualAdjustMonoSyn',false);

% GUI to manually curate cell classification. This is not neccesary at this point.
% It is more useful when you have multiple sessions
cell_metrics = CellExplorer('metrics',cell_metrics);

%% 3 - copy data to main session folder
% After you have run all this, you need to copy the Kilosort folder
% (actully only some files are neccesary), plus spikes.cellinfo.mat and cell_metrics.cellinfo.mat
% to the main session folder in the shared network drive

%% To only load spike times for futher analysis
[spikeT] = importSpikes('spikes',spikes,'UID',[]);

%% To load multiple sessions into CellExplorer
% After you have multiple sessions spike sorted and processed, you can open
% them all in cell explorer.

% assuming your current directory is (../project/animal/session), this will
% cd to project level (for example: A:\Data\GirardeauG)
cd('..\..')

% make this current directory into variable
data_path = pwd;

% look for all the cell_metrics.cellinfo.mat files
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);

% pull out basepaths and basenames
for i = 1:length(files)
    basepath{i} = files(i).folder;
    basename{i} = basenameFromBasepath(files(i).folder);
end

% write data frame of basepaths to project folder
% using this .csv, you will not need to again search for all the cell_metrics
% files, which can take some time if you have many sessions
df = table();
df.basepath = basepath';
df.basename = basename';
writetable(df,['sessions_',date,'.csv'])

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath,'basenames',basename);

% pull up gui to inspect all units in your project
cell_metrics = CellExplorer('metrics',cell_metrics);

% end

