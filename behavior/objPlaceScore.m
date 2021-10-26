function [objScore] = objPlaceScore(sessionsSeq,basepath,fs,timeTh)
% [objScore] = objPlaceScore(sessionsSeq,basepath)
% score memory performance in object location task, anotated online throug
% digital inputs

% sessionsSeq: sequence of recording sessions in temporal order
%    code:0 = sleep or other
%         1 = training
%         2 = test
%    default: [0 0 1 0 2] 
% fs default: 30000
% timeTh: [time training time test] to include in the analyis in seconds (discard time after it)
%    default: [300 300]
% 
% AntonioFR, 2020

if ~exist('basepath','var')
    basepath = pwd;
end
if ~exist('sessionsSeq','var')
    sessionsSeq = [0 0 1 0 2];
end   
if ~exist('fs','var')
    fs = 20000;
end    
if ~exist('timeTh','var')
    timeTh = [300 300];
end  

%% order individual ses
cd(basepath);
%basename = bz_BasenameFromBasepath(basepath);
%newdatpath = fullfile(basepath,[basename,'.dat']);
sespaths = [];
names2sort = [];
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
            end
        end
    end
end

[~,I] = sort(names2sort);
for idx = 1:length(I)
    nsespaths{idx} = sespaths{I(idx)};
end

% sanity check 
if length(sessionsSeq) ~= length(nsespaths)
    warning('inconsistent number of trials sessions');
end

%% get scores from each session
% objA=control; objB=displaced
cd([basepath '\' nsespaths{find(sessionsSeq == 1)}]);
data_training = getDigitalIn('all','fs', fs, 'offset', 0, 'periodLag', 1);
cd([basepath '\' nsespaths{find(sessionsSeq == 2)}]);
data_test = getDigitalIn('all','fs', fs, 'offset', 0, 'periodLag', 1);

cd(basepath);

%
if isempty(timeTh)
    object(1,:) = sum(data_training.dur{1,2}); % time in seconds spent with interacting object A during training
    object(2,:) = sum(data_training.dur{1,3}); % time in seconds spent with interacting object B during training
    object(3,:) = sum(data_test.dur{1,2}); %time in seconds spent with interacting object A during test
    object(4,:) = sum(data_test.dur{1,3}); %time in seconds spent with interacting object B during test
else
    times1 = find(data_training.timestampsOn{2}<timeTh(1));
    times2 = find(data_training.timestampsOn{3}<timeTh(1));
    object(1,:) = sum(data_training.dur{1,2}(times1)); 
    object(2,:) = sum(data_training.dur{1,3}(times2)); 
    
    times3 = find(data_test.timestampsOn{2}<timeTh(2));
    times4 = find(data_test.timestampsOn{3}<timeTh(2));
    object(3,:) = sum(data_test.dur{1,2}(times3)); 
    object(4,:) = sum(data_test.dur{1,3}(times4));     
end

object_A_training = object(1,:);
object_B_training = object(2,:);
object_A_test = object(3,:);
object_B_test = object(4,:);

% calculating discrimination index (novel-familiar / familiar + novel *100)
DI = (object_B_test - object_A_test) / (object_A_test + object_B_test) *100;
% calculating normalized object preference based on Ted Abel's 2014 paper
object_preference = object_B_test/(object_B_test + object_A_test)*100 - object_B_training/(object_B_training + object_A_training)*100;

%collect variables to save

objScore.object_training_time(:,1) = object_A_training;
objScore.object_training_time(:,2) = object_B_training;
objScore.object_test_time(:,1) = object_A_test;
objScore.object_test_time (:,2) = object_B_test;
objScore.object_test_time (:,3) = object_B_test + object_A_test;
objScore.object_training_time (:,3) = object_B_training + object_A_training;
objScore.discrimination_index = DI;
objScore.object_preference = object_preference;

% save variables to basepath
save ([basepath '\' 'objScore.mat'], 'objScore');

% figure;
% bar(object,0.5);
% set(gca,'XTickLabel',({'object A - novel', 'object B'}); % always check if proper object is indicated as novel!!
% ylabel('exploration time (s)');

% figure;
%     subplot(1,3,1);
%     bar(objScore.object_training_time,'r'); ylabel('time (s)'); set(gca,'XTickLabel',{'control', 'displaced'});hold on;
%     plot(trials.objNo,'b');hold on;xlim([1 length(trials.start)]);
%     ylabel('time (s)');xlabel('trial #');legend({'rew obj';'non rew obj'});title('obj exploration');
%     subplot(1,3,2);
%     plot(trials.expDif,'g');hold on;xlim([1 length(trials.start)]);
%     ylabel('time (s)');xlabel('trial #');title('dif in obj exploration');
%     subplot(1,3,3);
%     plot(1:1:length(auxT{1}),trials.durationAll,'b'); hold on;
%     plot(1:1:length(auxT{1}),trials.durationObj,'r'); hold on;xlim([1 length(trials.start)]);
%     ylabel('time (s)');xlabel('trial #');legend({'whole trial';'until last obj'});title('trial duration');
%     [path,name] = fileparts(basePath); title(name);



end

