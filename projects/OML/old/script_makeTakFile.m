X = BatchReturn('OMLCheese.batch');

for i=1:size(X,1)
    basepath = X{i,1};
    MergePoints = getStruct(basepath,'MergePoints');
    notsleep = cellfun(@(X) ~any(strfind(lower(X),'sleep')),MergePoints.foldernames);
    cd(basepath);
    folders = MergePoints.foldernames(notsleep);
    for j=1:length(folders) 
        cd(fullfile(basepath,folders{j}));
        q = dir('*.tak');
        if isempty(q) | q(1).bytes==0
            dlmwrite([folders{j}(1:end-14) '.tak'],'');
            delete([folders{1}(1:end-14) '.tak'],'');
        end
    end
end

