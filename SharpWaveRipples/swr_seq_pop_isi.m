function swr_seq_pop_isi(varargin)
%swr_seq_pop_isi - calculate ISIs during events such as SWRs between specfific cell types 

% Ryan Harvey, 10/21

p = inputParser;
addParameter(p,'basepath',[],@ischar) % single basepath to run 1 session
addParameter(p,'df',[]) % data frame with df.basepath
addParameter(p,'binary_class_variable','deepSuperficial',@ischar) % variable in cell_metrics that has groups
addParameter(p,'grouping_names',{'Deep','Superficial'},@iscell) % group names associated with the above
% cell array with group categories as strings (you still must specify binary_class_variable & grouping_names)
addParameter(p,'custom_grouping',{},@iscell) 
addParameter(p,'restrict_to_brainregion','CA1',@ischar) % brain region to run on (empty to run all)
addParameter(p,'restrict_to_celltype','Pyramidal Cell',@ischar) % cell class to run on (empty to run all)
addParameter(p,'force_run',false,@islogical) % to overwrite results
addParameter(p,'savepath',[]) % path to save results
addParameter(p,'parallel',true,@islogical) % run over sessions in parallel
addParameter(p,'shuffle',false,@islogical) % jitter spike times to make null dist
addParameter(p,'unique_unit_num',3,@isint) % min number of unique units per group per ripple
addParameter(p,'ripple_duration_restrict',[-inf,inf],@isnumeric) % min max duration of ripple
addParameter(p,'n_ripple_restrict',100,@isint) % min number of ripples for analysis
addParameter(p,'restrict_ts',[-inf,inf],@isnumeric) % restrict anaysis to time interval

parse(p,varargin{:})
basepaths = p.Results.basepath;
df = p.Results.df;
parallel = p.Results.parallel;
params = p.Results; % store other params for use in main

% % if empty basepaths, make it from sessions.csv (specific to Ryan H.)
% if isempty(basepaths)
%     df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
%     basepaths = unique(df.basepath);
% end

% convert basepaths to cell array if already not
if ~iscell(basepaths)
    basepaths = {basepaths};
end

% iter through main analysis. parallel or standard
WaitMessage = parfor_wait(length(basepaths));
if parallel
    parfor i = 1:length(basepaths)
        main(basepaths{i},params)
        WaitMessage.Send;
    end
else
    for i = 1:length(basepaths)
        main(basepaths{i},params)
        WaitMessage.Send;
    end
end
WaitMessage.Destroy;
end

function main(basepath,params)

disp(basepath)

file_str = strsplit(basepath,filesep);
savepath = fullfile(params.savepath,[file_str{end-1},'_',file_str{end},'.csv']);
if exist(savepath,'file') && ~params.force_run
    return
end

basename = basenameFromBasepath(basepath);

% load needed data
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']));
load(fullfile(basepath,[basename '.ripples.events.mat']));

% add custom_grouping to cell_metrics
if ~isempty(params.custom_grouping)
    cell_metrics.(params.binary_class_variable) = params.custom_grouping;
end

% make sure there are at least unique_unit_num of a each cell category available
good_to_run = check_unit_counts_per_group(cell_metrics,...
    params.restrict_to_brainregion,params.restrict_to_celltype,...
    params.binary_class_variable,params.grouping_names,params.unique_unit_num);
if ~good_to_run
    return
end

% restrict ripples to epoch
ripples = eventIntervals(ripples,params.restrict_ts,true);

% restrict ripples by duration
ripples = restrict_ripples_by_duration(ripples,params.ripple_duration_restrict);

% get ripple spikes
ripSpk = load_spikes_and_get_ripSpk(basepath,ripples,...
    params.restrict_to_brainregion,params.restrict_to_celltype);

% no units in ripples...skip
if ~isfield(ripSpk,'UnitEventAbs')
    return
end

% restict by number of unique units
ripSpk = restrict_ripples_unique_units(ripSpk,cell_metrics,params.grouping_names,...
    params.binary_class_variable,params.unique_unit_num);

% pass if fewer than allowed ripples
if length(ripSpk.EventRel) < params.n_ripple_restrict
    return
end

% get group isi distributions of each group and across groups
[A,B,AB] = calc_isi(ripSpk,cell_metrics,params.binary_class_variable,...
    params.grouping_names);

% get shuffled distributions
if params.shuffle
    [A_shuff,B_shuff,AB_shuff] = shuffle_isi(ripSpk,...
        cell_metrics,params.binary_class_variable,params.grouping_names);
else
    A_shuff = NaN;
    B_shuff = NaN;
    AB_shuff = NaN;
end
% save results as csv to savepath
save_results(A,B,AB,A_shuff,B_shuff,AB_shuff,params.grouping_names,savepath)
end

function good_to_run = check_unit_counts_per_group(cell_metrics,...
    restrict_to_brainregion,restrict_to_celltype,...
    binary_class_variable,grouping_names,unique_unit_num)

good_to_run = false;
if ~isempty(restrict_to_brainregion) && ~isempty(restrict_to_celltype)
    idx = contains(cell_metrics.brainRegion,restrict_to_brainregion) &...
        contains(cell_metrics.putativeCellType,restrict_to_celltype);
    
    if sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{1})) >= unique_unit_num &&...
            sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{2})) >= unique_unit_num
        good_to_run = true;
    end
elseif ~isempty(restrict_to_brainregion) % just restrict brain region
    idx = contains(cell_metrics.brainRegion,restrict_to_brainregion);
    
    if sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{1})) >= unique_unit_num &&...
            sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{2})) >= unique_unit_num
        good_to_run = true;
    end
elseif ~isempty(restrict_to_celltype) % just restrict cell type
    idx = contains(cell_metrics.putativeCellType,restrict_to_celltype);
    
    if sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{1})) >= unique_unit_num &&...
            sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{2})) >= unique_unit_num
        good_to_run = true;
    end
else % no restriction
    if sum(strcmp(cell_metrics.(binary_class_variable),grouping_names{1})) >= unique_unit_num &&...
            sum(strcmp(cell_metrics.(binary_class_variable),grouping_names{2})) >= unique_unit_num
        good_to_run = true;
    end
end
end

function ripples = restrict_ripples_by_duration(ripples,ripple_duration_restrict)

keep = ripples.duration >= ripple_duration_restrict(1) &...
    ripples.duration <= ripple_duration_restrict(2);

ripples.timestamps = ripples.timestamps(keep,:);
ripples.peaks = ripples.peaks(keep,:);
ripples.amplitude = ripples.amplitude(keep,:);
ripples.frequency = ripples.frequency(keep,:);
ripples.duration = ripples.duration(keep,:);
end

function ripSpk = load_spikes_and_get_ripSpk(basepath,ripples,...
    restrict_to_brainregion,restrict_to_celltype)

% if restricting brain and cell type
if ~isempty(restrict_to_brainregion) && ~isempty(restrict_to_celltype)
    spk = importSpikes('basepath',basepath,...
        'brainRegion',restrict_to_brainregion,...
        'cellType',restrict_to_celltype);
elseif ~isempty(restrict_to_brainregion) % just restrict brain region
    spk = importSpikes('basepath',basepath,...
        'brainRegion',restrict_to_brainregion);
elseif ~isempty(restrict_to_celltype) % just restrict cell type
    spk = importSpikes('basepath',basepath,...
        'cellType',restrict_to_celltype);
else % no restriction
    spk = importSpikes('basepath',basepath);
end

% make ripSpk struct with spike times per ripple
ripSpk = getRipSpikes('basepath',basepath,'events',ripples,'spikes',spk,...
    'saveMat',false);
end

function ripSpk = restrict_ripples_unique_units(ripSpk,cell_metrics,...
    grouping_names,binary_class_variable,unique_unit_num)

for i = 1:length(ripSpk.EventRel)
    [~,Locb] = ismember(ripSpk.EventRel{i}(2,:),cell_metrics.UID);
    keep(i) = sum(strcmp(cell_metrics.(binary_class_variable)(Locb), grouping_names{1})) >= unique_unit_num &&...
        sum(strcmp(cell_metrics.(binary_class_variable)(Locb), grouping_names{2})) >= unique_unit_num;
end

ripSpk.EventDuration = ripSpk.EventDuration(keep);
ripSpk.UnitEventAbs = ripSpk.UnitEventAbs(:,keep);
ripSpk.UnitEventRel = ripSpk.UnitEventRel(:,keep);
ripSpk.EventAbs = ripSpk.EventAbs(keep);
ripSpk.EventRel = ripSpk.EventRel(keep);
end

function [A,B,AB] = calc_isi(ripSpk,cell_metrics,binary_class_variable,...
    grouping_names)
% iter through each rip
A=[];   B=[];   AB=[];
for e = 1:numel(ripSpk.EventRel)
    for i = 1:size(ripSpk.EventRel{e},2)-1
        for j = 1:size(ripSpk.EventRel{e},2)-1
            [A,B,AB] = get_isi(ripSpk,...
                cell_metrics,...
                j,i,e,...
                A,B,AB,...
                binary_class_variable,...
                grouping_names);
        end
    end
end
end

function [A,B,AB] = get_isi(ripSpk,cell_metrics,j,i,e,A,B,AB,...
    binary_class_variable,grouping_names)
% get_isi: main function to calc isi on sub populations

if j ~= i && ripSpk.EventRel{e}(2,i) ~= ripSpk.EventRel{e}(2,j)
    if strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)})
        
        A=cat(1,A,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
        
    elseif strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)})
        
        B=cat(1,B,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
        
    elseif (strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)})) || ...
            (strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)}))
        
        AB=cat(1,AB,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
    end
end
end

function [A_shuff_all,B_shuff_all,AB_shuff_all] = shuffle_isi(ripSpk,...
    cell_metrics,binary_class_variable,grouping_names)

shuff_ripSpk = ripSpk;
% -0.01 to 0.01 seconds (10ms is ~ the longest ripple cycle)
min_s = -1/100;
max_s = 1/100;

A_shuff_all = [];
B_shuff_all = [];
AB_shuff_all = [];

% 10 shuffles is enough to de-couple cycle nesting of units
for shuff_i = 1:10
    for i = 1:length(ripSpk.EventRel)
        offset = (max_s-min_s).*rand(1,size(shuff_ripSpk.EventRel{i},2)) + min_s;
        shuff_ripSpk.EventRel{i}(1,:) = ripSpk.EventRel{i}(1,:) + offset;
    end
    [A_shuff,B_shuff,AB_shuff] = calc_isi(shuff_ripSpk,...
        cell_metrics,...
        binary_class_variable,...
        grouping_names);
    A_shuff_all = [A_shuff_all;A_shuff];
    B_shuff_all = [B_shuff_all;B_shuff];
    AB_shuff_all = [AB_shuff_all;AB_shuff];
end
end

function save_results(A,B,AB,A_shuff,B_shuff,AB_shuff,grouping_names,savepath)

% nan matrix
isi = nan(max([length(A),length(B),length(AB),...
    length(A_shuff),length(B_shuff),length(AB_shuff)]),6);

% fill each column of matrix
isi(1:numel(A),1) = A;
isi(1:numel(B),2) = B;
isi(1:numel(AB),3) = AB;
isi(1:numel(A_shuff),4) = A_shuff;
isi(1:numel(B_shuff),5) = B_shuff;
isi(1:numel(AB_shuff),6) = AB_shuff;

% column names
varnames = {grouping_names{1},...
    grouping_names{2},...
    [grouping_names{1},grouping_names{2}],...
    [grouping_names{1},'_shuff'],...
    [grouping_names{2},'_shuff'],...
    [grouping_names{1},grouping_names{2},'_shuff']};

% make table
df = table(isi(:,1),isi(:,2),isi(:,3),isi(:,4),isi(:,5),isi(:,6),...
    'VariableNames',...
    varnames);

writetable(df,savepath)
end