function mkSpks(savePath)
basepath = pwd;
basename = basenameFromBasepath(pwd);
if nargin <1
    savePath = basepath;
end
if exist([basename '.cell_metrics.cellinfo.mat'])
    load([basename '.cell_metrics.cellinfo.mat']);
else
    warning('Cell Metrics does not exist, computing');
    load([basename '.session.mat']);
    cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false); close all
end
if ~isfield(cell_metrics,'brainRegion')
    warning('Cell Metrics not fully calculated, recalculating')
    load([basename '.session.mat']);
    cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false); close all
end

useReg = unique(cell_metrics.brainRegion);
for i = 1:length(useReg)
    spikes = [];
    tempReg = []; tempReg = convertCharsToStrings(useReg{i});
    spikes = importSpikes('cellType', "Pyramidal Cell", 'brainRegion', tempReg);
    save([savePath '\' basename '.' useReg{i} 'pyr.cellinfo.mat'], 'spikes');
%     spikes = [];
%     spikes = importSpikes('cellType', ["Narrow Interneuron"; "Wide Interneuron"], 'brainRegion', tempReg);
%     save([savePath '\' basename '.' useReg{i} 'inter.cellinfo.mat'], 'spikes');
end
spikes = [];
spikes = importSpikes('cellType', "Pyramidal Cell");
save([savePath '\' basename '.allpyr.cellinfo.mat'], 'spikes');
spikes = [];
spikes = importSpikes('cellType', ["Narrow Interneuron"; "Wide Interneuron"]);
save([savePath '\' basename '.allint.cellinfo.mat'], 'spikes');
end