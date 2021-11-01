function [regID,cellID] = getSubs(cell_metrics)

regID = ones(length(cell_metrics.brainRegion),1);
check = ["CA1" "CA2" "CA3" "CTX" "DG"];

for i = 1:length(check)
    for j = 1:length(cell_metrics.brainRegion)
        if contains(cell2mat(cell_metrics.brainRegion(j)),convertStringsToChars(check(i)))
            regID(j) = (i+1);
        end
    end
end

cellID = ones(length(cell_metrics.putativeCellType),1);
inds = ismember(1:length(cellID),cell_metrics.tags.P);
cellID(inds) = 2;
inds = ismember(1:length(cellID),cell_metrics.tags.N);
cellID(inds) = 3;
end