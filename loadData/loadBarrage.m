function HSE = loadBarrage(basepath)
home_path = pwd;
cd([basepath '\Barrage_Files']);
basename = basenameFromBasepath(basepath);
anim = animalFromBasepath(basepath);

disp(strcat('Loading barrages from '," ",anim, "/", basename));
load([basename '.HSE.mat'],'HSE');
% load([basepath '\' basename '.session.mat']);

% pre=[]; post=[];
% for i=1:size(session.epochs,2)
%    if isempty(pre)
%        if contains(session.epochs{1,i}.name, 'pre')||contains(session.epochs{1,i}.name, 'Pre')
%            pre=i;
%        end
%    end
%    if isempty(post)
%        if contains(session.epochs{1,i}.name, 'post')||contains(session.epochs{1,i}.name, 'Post')
%            post=i;
%        end
%    end
% end
% 
% if isempty(pre)
%     for i=1:size(session.epochs,2)
%         if contains(session.epochs{1,i}.name, 'sleep')||contains(session.epochs{1,i}.name, 'homecage')
%              pre=i;
%         end
%     end
% end

% if isfield(session.epochs{1,pre}, 'startTime')
%     if isnumeric(session.epochs{1,pre}.startTime)
%         preAdjust = session.epochs{1,pre}.startTime;
%     else
%         preAdjust = str2double(session.epochs{1,pre}.startTime);
%     end
% else
%     preAdjust = 0;
% end
% 
% if isfield(session.epochs{1,post}, 'startTime')
%     if isnumeric(session.epochs{1,post}.startTime)
%         postAdjust = session.epochs{1,post}.startTime;
%     else
%         postAdjust = str2double(session.epochs{1,post}.startTime);
%     end
% else
%     error('Bad post adjustment');
% end
% 
% if (~isempty(HSE.pre.timestamps))&&(~isempty(HSE.post.timestamps))
%     correct_HSE.timestamps = [HSE.pre.timestamps+preAdjust; HSE.post.timestamps+postAdjust];
%     correct_HSE.peaks = [HSE.pre.peaks+preAdjust; HSE.post.peaks+postAdjust];
%     correct_HSE.duration = correct_HSE.timestamps(:,2)-correct_HSE.timestamps(:,1);
% 
%     clear HSE
%     HSE = correct_HSE;
%     HSE.keep = (1:size(HSE.timestamps,1)); %for this version Aza has imposed keep parameters separately
% 
%     load(strcat(basepath,'\',basename,'.SleepState.states.mat'));
%     SleepState.ints.NREMstate(:,1)=SleepState.ints.NREMstate(:,1)-1; %buffer?
%     SleepState.ints.NREMstate(:,2)=SleepState.ints.NREMstate(:,2)+1;
% 
%     if ~isempty(SleepState.ints.NREMstate)
%         HSEnREM = eventIntervals(HSE,SleepState.ints.NREMstate,1);
%         if ~isempty(HSEnREM)
%             [~,HSE.NREM] = intersect(HSE.peaks(HSE.keep), HSEnREM.peaks);
%         else
%             HSE.NREM = [];
%         end
%     end
% else
%     clear HSE
%     HSE.timestamps = []; HSE.peaks = []; HSE.duration = []; HSE.keep = []; HSE.NREM = [];
% end

cd(home_path);
end