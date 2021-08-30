function PreprocessSpikes(session) 
% This is a wrapper to concantenate the basic spike related functions,
% meant to be run right after the manual clustering

%% 1- compute spike stuff usign Peter scrips
    f = dir('Kilosort*');
    spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);

%% 2 - compute cell metrics using Peter scripts   
    cell_metrics = ProcessCellMetrics('session',session);
    
    % GUI to manually curate cell classification
    cell_metrics = CellExplorer('metrics',cell_metrics);
    
%% 3 - Load spike times 
    [spikeT] = bz_ImportSpikes('spikes',spikes,'UID',[]);
    
    
end

