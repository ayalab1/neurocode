function [regID,modID,regKey,modKey] = getSubs(cell_metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Subtypes
% Sort by region and modulation type for ALL cells. This will return
% indices corresponding to the regions, as well as the corresponding
% identification for these indices.
% 
% Should eventually make this more flexible for wider application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
inds = ismember(1:length(modID),cell_metrics.tags.P);
modID(inds) = 2;
inds = ismember(1:length(modID),cell_metrics.tags.N);
modID(inds) = 3;

regKey(1,:) = ["Unknown" check];
regKey(2,:) = 1:(length(check)+1);
modKey(1,:) = ["Unknown" "P" "N"];
modKey(2,:) = 1:3;
end