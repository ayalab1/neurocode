function [timestamps, idx, dayID] = pullMergePointsDay(basepath)
%Pull appropriate MergePoints timestamps to separate by day, in case there
%are subsessions within day

%MergePoints should be already ordered by date/time so we're just finding
%splitting points here

% Outputs:
% timestamps: new compiled start/stop timestamps per day (Nx2 where N is
%             number of days
% idx: MergePoint idxs pulled for the day indicating the first and last
%             sessions of the day
% dayID: overarching basename day number based on folder information. This
% functionality is dependent upon having a structured animal/basename
% format where the basename includes "day#". 

basename = basenameFromBasepath(basepath);
load([basepath '\' basename '.MergePoints.events.mat']);
fun_existsdate = @(x) (strfind(ismember(x, '1234567890'), [0, ones(1, 6), 0])+1); %locate date information
dateStore = [];
for i = 1:size(MergePoints.foldernames,2)
    getIDX = fun_existsdate(MergePoints.foldernames{1,i});
    dateStore(i) = str2num(MergePoints.foldernames{1,i}(getIDX:(getIDX+5)));
end
timestamps = [];
idx = [];
getUnique = unique(dateStore);
for i = 1:length(getUnique)
    [~,dateIDX] = find(dateStore==getUnique(i));
    if length(dateIDX)==1
        timestamps = [timestamps; MergePoints.timestamps(dateIDX,:)];
        idx = [idx; dateIDX, dateIDX];
    else
        timestamps = [timestamps; MergePoints.timestamps(min(dateIDX),1), MergePoints.timestamps(max(dateIDX),2)];
        idx = [idx; min(dateIDX), max(dateIDX)];
    end
end

if nargout==3
    dayID = getDayNum(getUnique);
end
end

%%
function [dayIDs] = getDayNum(getUnique)
%helper function for pullMergePointsDay to identify the reference day
%number for each identified day
previousDir = pwd;
cd ..
allFolders = dir();
fun_existsday = @(x) (strfind(ismember(x, 'Day')|ismember(x, 'day'), [ones(1,3)])); %locate date information
fun_existsSubfolder = @(x) (strfind(ismember(x, '1234567890'), [0, ones(1, 6), 0, ones(1,6)])+1); %locate subfolder, NOT KILOSORT
fun_existsKiloSort = @(x) (strfind(ismember(x, 'Kilosort'), [ones(1,8)])); %locate date information
fun_existsNumber = @(x) (strfind(ismember(x, '1234567890'), [ones(1, 1)])); %locate day numbers in folder name
found = false; 
i = 0;
while ~found
    i=i+1; %control counter
    if i <= size(allFolders,1)
        if ~isempty(fun_existsday(allFolders(i).name))
            useIDX = i;
            cd([allFolders(useIDX).folder '\' allFolders(useIDX).name]);
            subDir = dir();
            for j = 1:size(subDir,1)
                if ~isempty(fun_existsSubfolder(subDir(j).name))&&isempty(fun_existsKiloSort(subDir(j).name))
                    startIDX = fun_existsSubfolder(subDir(j).name);
                    useIDXsub = j;
                    found = true;
                end
            end
        end
    else
        error('No day folder found, check');
    end
end

useDate = subDir(useIDXsub).name(startIDX:(startIDX+5));
refDate = datetime(str2num(['20' useDate(1:2)]), str2num(useDate(3:4)), str2num(useDate(5:6)));
dayDiff = [];
for i = 1:size(getUnique,2)
    tempStr = num2str(getUnique(i));
    tempDate = datetime(str2num(['20' tempStr(1:2)]), str2num(tempStr(3:4)), str2num(tempStr(5:6)));
    dayDiff(i) = daysdif(refDate,tempDate);
end

useDay = allFolders(useIDX).name;
getNums = fun_existsNumber(useDay);
flagging = diff(getNums)>1;
if size(flagging,1)==1
    refDay = str2num(useDay(getNums));
else
    error('Not implemented');
end

dayIDs = dayDiff+refDay;
cd(previousDir);
end