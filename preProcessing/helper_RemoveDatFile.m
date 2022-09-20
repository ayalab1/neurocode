%  helper_RemoveDatFile- script that checks if the .dat file in the current directory matches the size of the individual subsession .dat files and removes it

clear basepath basename bytes q MergePoints folders

f = dir('Kilosort*');
if isempty(f), display(['No Kilosort folder! Aborting...']), return; end

basepath = pwd;
basename = basenameFromBasepath(basepath);
MergePoints = getStruct(basepath,'MergePoints');
folders = MergePoints.foldernames';
for i=1:length(folders),
    folder = folders{i};
    q = dir(fullfile(folder,'amplifier.dat'));
    try
        bytes(i,1) = q.bytes;
    catch
        bytes(i,1) = 0;
    end
end

q = dir(fullfile(basepath,[basename '.dat']));
if isempty(q), display([basename '.dat is already removed!']);
elseif sum(bytes)==q.bytes,
    delete(fullfile(basepath,[basename '.dat']));
    display([basename '.dat deleted!']);
else
    display(['Subsession dat files don''t match ' basename '.dat. Aborting..']);
end