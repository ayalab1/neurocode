%  helper_RemoveDatFile- script that checks if the .dat file in the current directory matches the size of the individual subsession .dat files and removes it

% EXAMPLE: 
% % To remove all the concatenated .dat files in the OJR project, I ran the following code:
% projectFolder = 'Y:\OJRproject'; % this is the only line you need to change to run the script in all the folders of your project
% animalFoldersStruct = dir(projectFolder); % find all the subfolders in the project folder
% animalFolders = struct2cell(animalFoldersStruct); 
% animalFolders = animalFolders(1,Unfind(3:length(animalFoldersStruct)) & cat(1,animalFoldersStruct.isdir))'; % convert to a list of folders
% for animal = 1:length(animalFolders) % for each animal
%     sessionFoldersStruct = dir(fullfile(projectFolder,animalFolders{animal})); % list the folders containing recording sessions
%     sessionFolders = struct2cell(sessionFoldersStruct); sessionFolders = sessionFolders(1,Unfind(3:length(sessionFoldersStruct)) & cat(1,sessionFoldersStruct.isdir))'; % convert to a list of folders
%     for pathID = 1:length(sessionFolders) % for each session
%         cd(fullfile(projectFolder, animalFolders{animal}, sessionFolders{pathID})); % cd to the session
%         try
%             helper_RemoveDatFile % run script to perform checks and remove concatenated .dat file
%         catch
%             warning(['something didn''t work with ' fullfile(projectFolder, animalFolders{animal}, sessionFolders{pathID})]);
%         end
%     end
% end

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