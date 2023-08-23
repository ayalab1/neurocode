function [objScore] = preferenceScore(basepath, timeTh, runDigIn, digChans)

if ~exist('basepath','var')
    basepath = pwd;
end
if ~exist('timeTh','var')
    timeTh = [300 300];
end  
if ~exist('runDigIn','var')
    runDigIn=false;
end
if ~exist('digChans','var')
    digChans = [2,3]; %channels to compare (A,B)
end

%%
if runDigIn
    cd(basepath);
    processFolder;
end

%% order individual ses
cd(basepath);
sespaths = [];
names2sort = [];
dates2sort = [];
nsespaths = [];
d = dir;
fidx = 0;

for idx = 3:length(d)
    if d(idx).isdir
        if (numel(d(idx).name) > 13)
            if(numel(num2str(str2num(d(idx).name(end-5:end))))>=5 && numel(num2str(str2num(d(idx).name(end-12:end-7))))==6) %detecting intan recordings
                fidx = fidx+1;
                sespaths{fidx} = fullfile(d(idx).name);
                names2sort(fidx) = str2num(d(idx).name(end-5:end));
                dates2sort(fidx) = str2num(d(idx).name(end-12:end-7));
            end
        end
    end
end
per_day = unique(dates2sort);
per_day = sort(per_day, 'ascend');
count = 1;
for i = 1:length(per_day)
    use_idx = (dates2sort==per_day(i));
    [~,sortBy] = sort(names2sort(use_idx==1));
    for j = 1:sum(use_idx)
        path_idx = find(use_idx==1);
        nsespaths{count} = sespaths{path_idx(sortBy==j)};
        count = count+1;
    end
end

%% Pull data per folder
trainCt= 1;
for i = 1:length(nsespaths)
    cd([basepath '\' nsespaths{i}]);
    try
        load('digitalIn.events.mat');
    catch
        warning('no Digital In detected, running processFolder');
        cd(basepath);
        processFolder;
        cd([basepath '\' nsespaths{i}]);
        load('digitalIn.events.mat');
    end
    if contains(nsespaths{i}, 'rain') %train or Train
       if isempty(timeTh)
           train_time(trainCt,1) = sum(digitalIn.dur{1,digChans(1)});
           train_time(trainCt,2) = sum(digitalIn.dur{1,digChans(2)});
       else
          clicker_trainA = find(digitalIn.timestampsOn{2}<timeTh(1));
          clicker_trainB = find(digitalIn.timestampsOn{3}<timeTh(1));
          train_time(trainCt,1) = sum(digitalIn.dur{1,digChans(1)}(clicker_trainA));
          train_time(trainCt,2) = sum(digitalIn.dur{1,digChans(2)}(clicker_trainB));
       end
       trainCt = trainCt+1; clear digitalIn;
    elseif contains(nsespaths{i}, 'test')
        if isempty(timeTh)
           test_time(1) = sum(digitalIn.dur{1,digChans(1)});
           test_time(2) = sum(digitalIn.dur{1,digChans(2)});
       else
          clicker_testA = find(digitalIn.timestampsOn{digChans(1)}<timeTh(1));
          clicker_testB = find(digitalIn.timestampsOn{digChans(2)}<timeTh(1));
          test_time(1) = sum(digitalIn.dur{1,digChans(1)}(clicker_testA));
          test_time(2) = sum(digitalIn.dur{1,digChans(2)}(clicker_testB));
       end
       trainCt = trainCt+1; clear digitalIn;
    end
end

DI = ((test_time(2)-test_time(1))/(test_time(2)+test_time(1)))*100;
obj_pref = (test_time(2)/(test_time(2)+test_time(1))*100) - (sum(train_time(:,2))/sum(train_time,'all')*100);

objScore.object_training_time = train_time;
objScore.object_test_time = test_time;
objScore.discrimination_index = DI;
objScore.object_preference = obj_pref;

save([basepath '\' 'objScore.mat'], 'objScore');
cd(basepath);
end