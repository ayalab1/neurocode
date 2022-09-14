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
writetable(df,['sessions_',date,'.csv']);

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath,'basenames',basename);

% pull up gui to inspect all units in your project
cell_metrics = CellExplorer('metrics',cell_metrics);