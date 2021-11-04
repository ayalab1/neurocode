function [regID,modID,regKey,modKey] = getSubs(cell_metrics, regCheck, tagCheck)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Subtypes
% Sort by region and modulation type for ALL cells. This will return
% indices corresponding to the regions, as well as the corresponding
% identification for these indices.
% 
% Should eventually make this more flexible for wider application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    tagCheck = ["P" "N"];
    if nargin < 2
        regCheck = ["CA1" "CA2" "CA3" "CTX" "DG"];
    end
end

regID = ones(length(cell_metrics.brainRegion),1);
check = ["CA1" "CA2" "CA3" "CTX" "DG"];

for i = 1:length(check)
    for j = 1:length(cell_metrics.brainRegion)
        if contains(cell2mat(cell_metrics.brainRegion(j)),convertStringsToChars(check(i)))
            regID(j) = (i+1);
        end
    end
end

modID = ones(length(cell_metrics.putativeCellType),1);
for i = 1:length(tagCheck)
    inds = ismember(1:length(modID),eval(strcat('cell_metrics.tags.',tagCheck(i))));
    modID(inds) = i+1;
end

regKey(1,:) = ["Unknown" regCheck];
regKey(2,:) = 1:(length(regCheck)+1);
modKey(1,:) = ["Unknown" tagCheck];
modKey(2,:) = 1:(length(tagCheck)+1);
end