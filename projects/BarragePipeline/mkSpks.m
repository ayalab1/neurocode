function mkSpks(savePath)
basepath = pwd;
basename = basenameFromBasepath(pwd);
if nargin <1
    savePath = basepath;
end
load([basename '.cell_metrics.cellinfo.mat']);
useReg = unique(cell_metrics.brainRegion);
for i = 1:length(useReg)
    spikes = [];
    tempReg = []; tempReg = convertCharsToStrings(useReg{i});
    spikes = importSpikes('cellType', "Pyramidal Cell", 'brainRegion', tempReg);
    save([savePath '\' basename '.' useReg{i} '.cellinfo.mat'], 'spikes');
end
spikes = [];
spikes = importSpikes('cellType', "Pyramidal Cell");
save([savePath '\' basename '.allpyr.cellinfo.mat'], 'spikes');
end