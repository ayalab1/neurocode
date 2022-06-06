%% Put info needed into a hmmLinearized structure

% load csv back to matlab
basepath = 'M:\Data\Can\OML22\day18';
basename = basenameFromBasepath(basepath);

csv = 'linearizedPosition.csv';
opts = detectImportOptions(fullfile(basepath,csv),'NumHeaderLines',1);
df_linearizedPos = readtable(fullfile(basepath,csv),opts);
df_linearizedPos = df_linearizedPos{:,:}; % make it into a matrix

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

hmmLinearized.linearized = df_linearizedPos(:,2);
hmmLinearized.projection.x = df_linearizedPos(:,4);
hmmLinearized.projection.y = df_linearizedPos(:,5);
hmmLinearized.timestamps = behavior.timestamps;
hmmLinearized.segment = df_linearizedPos(:,3);
hmmLinearized.nodes = node; % orginal 2d coord of nodes from clicking on matlab figure
hmmLinearized.linearizedBoundaries = df_boundCoord(:,2:3); % linearized boundaries of segments
hmmLinearized.edgeOrder = df_edgeOrder(:,2:3); %begin and end node number of each segment; 0 based from python

save(fullfile(basepath,[basename, '_hmmLinearized.mat']),'hmmLinearized');

%% Load the files we will be working with

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath); 

load(fullfile(basepath,[basename '_hmmLinearized_2.mat']));
behavior = getStruct(basepath,'animal.behavior');

plotFigures = true; % set to false to execute this script without plotting any figures

%% Merge "left" and "right" linearized positions (stretch "right" into "left")

l = hmmLinearized.linearized;
t = hmmLinearized.timestamps(:);
edges = hmmLinearized.linearizedBoundaries;
% keep the positions from the first 3 edges (of hmmLinearized.linearizedBoundaries(1:2,:))
% keep the positions from leftbound trials (edges 3:13)
edgesLeft = edges(3:13,:);
edgesRight = edges(14:24,:);
% interpolate positions from rightbound trials as if they were left bound trials:
% fit edges 14:24 to edges 3:13
right = l>edges(14,1);
l(right) = interp1([edgesRight(:,1); edgesRight(end,2)], [edgesLeft(:,1); edgesLeft(end,2)], l(right));

%% Define trials

laps = unwrap(l/max(l)*2*pi)/(2*pi);
boundaries = [0:floor(max(laps))-1]'; boundaries(:,2) = boundaries+1; % Note: we use floor(max(laps)) because the the last trial is the last finished trial (unfinished trials are not counted)
nTrials = size(boundaries,1);
trials = nan(nTrials,2);
% in reality, sometimes the rat would finish a trial and take some time before starting on the next trial
% he spends this time around the start port, i.e. changing the trial back and forth from (i-1) to (i)
% To avoid this time of non-performance (which can, on some occasions last for many minutes) as a part of any trial,
% one may define trial (i) to begin at a timepoint after which the rat does not go back to the positions of the previos
% trial (but advances forward):
for i=1:nTrials
    start = t(1+find(laps<boundaries(i,1),1,'last')); % just after the last moment the rat traverses into the current lap for the last time 
    stop = t(find(laps>boundaries(i,2),1,'first')-1); % just before the first moment the rat traverses into the next lap for the first time
    trials(i,:) = [start stop];
end

if plotFigures,
    figure; plot(t,laps);
    PlotHVLines(0:nTrials,'h','k');
    for i=1:nTrials
        PlotIntervals(trials(i,:),'color','r','ylim',boundaries(i,:));
    end
    ylabel('unwrapped position (lap by lap');
    legend('trials');
    xlabel('time (s)');
    ylim([min(laps) max(laps)]);
end

%%  assign trial ID (left vs right)

segmentID = zeros(size(t));
% assign "left" vs "right" identity to the segments. While all segments within
% "edgesLeft" are technically left, the first 3 segments (close to the choice point)
% and the last segment (close to the start box) are often assigned to the wrong direction
% The decision will therefore be made based on the 7 segments on the middle of the "left"
% and "right" trajectories (edgesLeft(4:end-1,:))
segmentID(InIntervals(hmmLinearized.linearized,edgesLeft(4:end-1,:))) = -1; % left
segmentID(InIntervals(hmmLinearized.linearized,edgesRight(4:end-1,:))) = 1; % right

[in,w] = InIntervals(t,trials); in = in & segmentID~=0;
direction = Accumulate(w(in),segmentID(in),'mode','mean');

i=0;
%% Human check: execute this for each trial

if plotFigures,
    i = i+1; % go on to the next trial
    this = segmentID;
    clf
    in = InIntervals(t,trials(i,:));
    ax = subplot(1,2,1);
    plot(behavior.position.x(1:10:end),behavior.position.y(1:10:end),'k.','markersize',1); hold all; % plot all the points as background
    scatter(behavior.position.x(in),behavior.position.y(in),30,t(in)-min(t(in)),'filled');
    set(get(colorbar,'YLabel'),'String','time in trial (s)','fontsize',12);
    colormap(gca,parula);
    title(['Trial ' num2str(i)]);

    ay = subplot(1,2,2);
    plot(behavior.position.x(1:10:end),behavior.position.y(1:10:end),'k.','markersize',1); hold all; % plot all the points as background
    scatter(hmmLinearized.projection.x(in)+randn(sum(in),1)*5,hmmLinearized.projection.y(in)+5*randn(sum(in),1),30,this(in),'filled');
    clim([-1 1]); clabel('-1 = left; 0 = ignored; 1 = right');
    colormap(gca,[155 199 255; 212 152 255; 255 145 145]/255);
    colormap(gca,[72 151 252; 233 203 255; 255 56 56]/255);
    set(colorbar,'ytick',[-1 0 1],'yticklabel',{'Left','Ignored','Right'},'fontsize',12);
    % EquateScales
    if direction(i)>0
        title('right')
    else
        title('left');
    end

    drawnow
end

%% Convert linearized positions from pixels to cm!!!

maze_sizes = 360; % length of one loop on figure8 maze is 360 cm (8 feet)
pos_range = max(l)-min(l);
cmPerPixel = maze_sizes/pos_range;
l = l*cmPerPixel;

%% Save trials and linearized positions

behavior.trials = trials;
behavior.trialID = direction; % -1 = left; 1 = right;
behavior.position.linearized = l';

save(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');















