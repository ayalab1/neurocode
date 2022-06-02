%% Put info needed into a hmmLinearized structure

% load csv back to matlab
basepath = 'M:\Data\Can\OML22\day18';
basename = basenameFromBasepath(basepath);

csv = 'linearizedPosition.csv';
opts = detectImportOptions(fullfile(basepath,csv),'NumHeaderLines',1);
df_linPos = readtable(fullfile(basepath,csv),opts);
df_linPos = df_linPos{:,:}; % make it into a matrix

csv = 'edgeOrder.csv';
opts = detectImportOptions(fullfile(basepath,csv),'NumHeaderLines',1);
df_edgeOrder = readtable(fullfile(basepath,csv),opts);
df_edgeOrder = df_edgeOrder{:,:}; % make it into a matrix

csv = 'bound_coord.csv';
opts = detectImportOptions(fullfile(basepath,csv),'NumHeaderLines',1);
df_boundCoord = readtable(fullfile(basepath,csv),opts);
df_boundCoord = df_boundCoord{:,:}; % make it into a matrix

% load node coordinates
load('node.mat')

hmmLinearized.linearized = df_linPos(:,2);
hmmLinearized.projection.x = df_linPos(:,4);
hmmLinearized.projection.y = df_linPos(:,5);
hmmLinearized.timestamps = behavior.timestamps;
hmmLinearized.segment = df_linPos(:,3);
hmmLinearized.nodes = node; % orginal 2d coord of nodes from clicking on matlab figure
hmmLinearized.linBoundaries = df_boundCoord(:,2:3); % linearized boundaries of segments
hmmLinearized.edgeOrder = df_edgeOrder(:,2:3); %begin and end node number of each segment; 0 based from python

save(fullfile(basepath,[basename, '_hmmLinearized.mat']),'hmmLinearized');
