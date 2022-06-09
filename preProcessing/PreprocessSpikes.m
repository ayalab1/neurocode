% function preprocessSpikes(session)
% This is a wrapper to concantenate the basic spike related functions,
% meant to be run right after the manual clustering

%% 0 - Spike sorting
% Automatic spike sorting should be run with Kilosort (using preprocessSession)
% After that, copy yo your local SSD the Kilosort folder, .dat and .xml files
% Do manual curation using Phy2

%% 1- extract spike times and waveforms for sorted clusters
%if you intend to concatenate multiple kilosort runs into one spikes structure
ifMultiKiloSort = 0; 
session = sessionTemplate(pwd,'showGUI',false);
f = dir('Kilosort*');
% Make sure there is only one KiloSort folder before running, unless you
% needed to spike sort probes separately (ifMultiKiloSort=1).

% The hippocampal KiloSort folder should be listed first in the session 
% folder for organization, but this is not necessary for running

% ADDITIONALLY: If combining multiple KiloSort files, MAKE SURE ALL
% SPIKE GROUPS IN NEUROSCOPE ARE REACTIVATED (in yellow tab in neuroscope,
% groups need to be assigned, no '?'). 

% check if spikes.cellinfo has already been created
pre_exist_spike_files = dir('*spikes*.cellinfo.mat');
if (size(f,1) > 1) && ifMultiKiloSort
    for i = 1:size(f,1)
        session.spikeSorting{1,1}.relativePath = f(i).name;
        
        % if all the spike files exist, load instead of making from scratch
        if length(pre_exist_spike_files) == length(f)
            temp_spikes = load(pre_exist_spike_files(i).name);
            tempSpikes{i} = temp_spikes.spikes;
        else
            tempSpikes{i} = loadSpikes('session',session,...
                'clusteringpath',[f(i).folder filesep f(i).name],...
                'savemat', false);
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
    basename = basenameFromBasepath(pwd);
    save(strcat(basename,'.spikes.cellinfo.mat'),'spikes');
    clear tempSpikes i j tempSpindices2 prevUID base comp
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

