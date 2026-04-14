function useReg = pullSpikes(basepath)
% Create spike structures per region and cell type
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% basepath:     Full path where session is located. Default: pwd

%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% useReg:       Cell containing strings of the regions present in cell
%               metrics, used for loading in some BARR-related analyses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Handle inputs
if nargin<1
    warning('No basepath provided, defaulting to current directory');
    basepath = pwd;
end

original = pwd;
cd(basepath);
basename = basenameFromBasepath(basepath);

cd(basepath);
%% Load and prepare cell_metrics
if exist([basepath filesep basename '.cell_metrics.cellinfo.mat'])
    load([basepath filesep basename '.cell_metrics.cellinfo.mat']);
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
    save([basepath filesep basename '.cell_metrics.cellinfo.mat'], 'cell_metrics');
end

%% Pull regions
tempUseReg = unique(cell_metrics.brainRegion);
useReg = []; tuC = 1; y1 = 0; y2 = 0; y3 = 0;
for i = 1:length(tempUseReg)
    if contains(convertCharsToStrings(tempUseReg{i}),"CA1") %in case sublayers labeled, collapse all into CA1
        if ~y1
            useReg{tuC} = "CA1";
            tuC = tuC+1; y1=1;
        end
    elseif contains(convertCharsToStrings(tempUseReg{i}),"CA2")
        if (~y2)
            useReg{tuC} = "CA2";
            tuC = tuC+1; y2=1;
        end
    elseif contains(convertCharsToStrings(tempUseReg{i}),"CA3")
        if (~y3)
            useReg{tuC} = "CA3";
            tuC = tuC+1; y3=1;
        end
    else
        useReg{tuC} = convertCharsToStrings(tempUseReg{i});
        tuC = tuC+1;
    end
end

cd(original);
end