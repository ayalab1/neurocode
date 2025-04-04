function pullSpikes(basepath, savePath)
% Create spike structures per region and cell type
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% basepath:     Full path where session is located. Default: pwd
% savePath:     Location for spike structures to be saved. For example,
%               '[basepath '\Barrage_Files']'. Default: pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Handle inputs
if nargin < 1
    basepath = pwd;
    savePath = basepath;
elseif nargin <2
    savePath = basepath;
end
original = pwd;
cd(basepath);
basename = basenameFromBasepath(pwd);

%% Load and prepare cell_metrics
if exist([basepath '\' basename '.cell_metrics.cellinfo.mat'])
    load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
else
    warning('cell_metrics does not exist, computing');
    load([basename '.session.mat']);
    cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false); close all
    clear session
end

if ~isfield(cell_metrics,'brainRegion')
    error('cell_metrics incomplete, recalculate and rerun');
end

if ~isfield(cell_metrics,'tags')
    cell_metrics.tags.Bad = [];
    save([basepath '\' basename '.cell_metrics.cellinfo.mat']);
end

%% Pull regions
tempUseReg = unique(cell_metrics.brainRegion);
useReg = []; tuC = 1; y1 = 0; y2 = 0; y3=0;
for i = 1:length(tempUseReg)
    if contains(convertCharsToStrings(tempUseReg{i}),"CA1")
        if ~y1
            useReg{tuC} = "CA1";
            tuC = tuC+1; y1 = 1;
        end
    elseif contains(convertCharsToStrings(tempUseReg{i}),"CA2")
        if (~y2)
            useReg{tuC} = "CA2";
            tuC = tuC+1; y2 = 1;
        end
    elseif contains(convertCharsToStrings(tempUseReg{i}),"CA3")
        if (~y3)
            useReg{tuC} = "CA3";
            tuC = tuC+1; y3 = 1;
        end
    else
        useReg{tuC} = convertCharsToStrings(tempUseReg{i});
        tuC = tuC+1;
    end
end

%% Pull pyramidal per region
for i = 1:length(useReg)
    spikes = [];
    spikes = importSpikes('cellType', "Pyramidal Cell", 'brainRegion', useReg{i});
    countCheck = 0; UIDcheck = [];
    for m = 1:length(cell_metrics.brainRegion)
        if contains(cell_metrics.brainRegion{m}, useReg{i})
            if contains(cell_metrics.putativeCellType{m},'Pyr')
                if isfield(cell_metrics.tags,'Bad')
                    if ~ismember(cell_metrics.UID(m),cell_metrics.tags.Bad)
                        countCheck = countCheck+1;
                        UIDcheck = [UIDcheck cell_metrics.UID(m)];
                    end
                else
                    countCheck = countCheck+1;
                    UIDcheck = [UIDcheck cell_metrics.UID(m)];
                end
            end
        end
    end
    if sum(UIDcheck==spikes.UID)==length(spikes.times)
        save(strcat(savePath, '\', basename, '.', useReg{i}, 'pyr.cellinfo.mat'), 'spikes');
    else
        error('Issue with importing');
    end
    
    %% Pull interneurons per region
    spikes = [];
    spikes = importSpikes('cellType', "Int", 'brainRegion', useReg{i});
    countCheck = 0; UIDcheck = [];
    for m = 1:length(cell_metrics.brainRegion)
        if contains(cell_metrics.brainRegion{m}, useReg{i})
            if contains(cell_metrics.putativeCellType{m},'Int')
                if isfield(cell_metrics.tags,'Bad')
                    if ~ismember(cell_metrics.UID(m),cell_metrics.tags.Bad)
                        countCheck = countCheck+1;
                        UIDcheck = [UIDcheck cell_metrics.UID(m)];
                    end
                else
                    countCheck = countCheck+1;
                    UIDcheck = [UIDcheck cell_metrics.UID(m)];
                end
            end
        end
    end
    if sum(UIDcheck==spikes.UID)==length(spikes.times)
        save(strcat(savePath, '\', basename, '.', useReg{i}, 'int.cellinfo.mat'), 'spikes');
    else
        error('Issue with importing');
    end
end

%% Pull all pyr or int cells
spikes = [];
spikes = importSpikes('cellType', "Pyramidal Cell");
save([savePath '\' basename '.allpyr.cellinfo.mat'], 'spikes');
spikes = [];
spikes = importSpikes('cellType', ["Narrow Interneuron"; "Wide Interneuron"]);
save([savePath '\' basename '.allint.cellinfo.mat'], 'spikes');

cd(original);
end