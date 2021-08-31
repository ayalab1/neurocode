function preprocessSpikes(session) 
% This is a wrapper to concantenate the basic spike related functions,
% meant to be run right after the manual clustering

%% 0 - Spike sorting
  % Automatic spike sorting should be run with Kilosort (using preprocessSession)
  % After that, copy yo your local SSD the Kilosort folder, .dat and .xml files
  % Do manual curation using Phy2
  
%% 1- extract spike times and waveforms for sorted clusters
    session = sessionTemplate(pwd,'showGUI',false);
    f = dir('Kilosort*');
    spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);

%% 2 - compute basic cell metrics 
    cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
   
    % GUI to manually curate cell classification. This is not neccesary at this point. 
    % It is more useful when you have multiple sessions
    cell_metrics = CellExplorer('metrics',cell_metrics);
    
%% 3 - copy data to main session folder
   % After you have run all this, you need to copy the Kilosort folder
   % (actully only some files are neccesary), plus spikes.cellinfo.mat and cell_metrics.cellinfo.mat
   % to the main session folder in the shared network drive
    
%% To only load spike times for futher analysis
    [spikeT] = importSpikes('spikes',spikes,'UID',[]);
    
    
end

