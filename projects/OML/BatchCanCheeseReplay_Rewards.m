function [  scores,averages,outputsCurrent,outputsPrev,pre,post,durations,distanceThreshold,  conditionType,taskType,basepath] = BatchCanCheeseReplay_Rewards(basepath,conditionType,taskType)

%%
if conditionType>1 || taskType>1, return; end

options = [5 10 20 25 30 40 50 Inf];
optionIndex = 2; % use only this index

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
cd(basepath);
behaviorFileChanged = false;

disp([datestr(clock) ': Starting session ' basepath '...']);

try
    behavior = getStruct(basepath,'animal.behavior');
catch
    disp('No behavior file!')
    error
end
try
    trials = behavior.trials;
    preProbeTrials = behavior.preProbeTrials;
    if isempty(trials), error; end
catch
    load('cheesbTrialsDayOffset.mat')
    trials = [cheesbTrialsDay.trials.start cheesbTrialsDay.trials.end];
    behavior.trials = trials;
    preProbeTrials = [cheesbTrialsDay.preProbe.start cheesbTrialsDay.preProbe.end];
    behavior.preProbeTrials = preProbeTrials;
    behavior.postProbeTrials = [cheesbTrialsDay.postProbe.start cheesbTrialsDay.postProbe.end];
    
end

% trials are sometimes wrong. Take the whole epoch
MergePoints = getStruct(basepath,'MergePoints');
trials = ConsolidateIntervals(MergePoints.timestamps(CountInIntervals(trials(:,1),MergePoints.timestamps)>0,:),'epsilon',1);
% pre/post probe _epoch_ is sometimes also wrong (e.g. OML5day6). Take all running before/after trials
preProbeTrials = [0 trials(1)];
postProbeTrials = [trials(end) Inf];

try % check if linearized positions are saved
    l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
    l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
    boundaries = behavior.position.currentRewards.boundaries;
    boundaries = behavior.position.previousRewards.boundaries;
   
catch 
    behaviorFileChanged = true;
    % compute linearized positions relative to current or past rewards
    % reward configurations:
     
    rewardsConfig{6} = [5,8; 9,4; 11, 12];
    rewardsConfig{7} = [5 ,4; 11,8; 8 ,13];
    rewardsConfig{8} = [5, 13; 8, 8; 12, 4];
    rewardsConfig{9} = [3, 5; 10, 5; 9, 13];
    rewardsConfig{10} = [5, 11; 13, 9; 8, 3];
    rewardsConfig{11} = [5, 4; 8, 8; 14, 4];
    rewardsConfig{12} = [4, 8; 12, 7; 10, 15];
    rewardsConfig{13} = [3, 3; 8, 10; 14, 10];
    rewardsConfig{14} = [6, 8; 11, 5; 5, 14];
    rewardsConfig{15} = [7, 2; 13, 7; 9, 12];
    rewardsConfig{16} = [4, 11; 13, 13; 8, 6];
    rewardsConfig{17} = [4, 4; 12, 8; 7, 14];
    rewardsConfig{18} = [4, 12; 8, 8; 11, 3];
    % each OML animal is a cell. The cell contains a lookup table of [day,reward configuration]
    rewardOrder{1} = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 11; 12 12; 13 13; 14 14; 15 15; 16 16; 17 17; 18 18; 19 nan];
    rewardOrder{3} = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 11; 12 12; 13 13; 14 14; 15 15; 16 16; 17 17; 18 18; 19 nan];
    rewardOrder{5} = [1 6; 2 7; 3 8; 4 9; 5 10; 6 11; 7 12; 8 13; 9 14; 10 15; 11 nan; 12 16; 13 17; 14 nan];
    rewardOrder{7} = [1 6; 2 7; 3 8; 4 9; 5 11; 6 11; 7 12; 8 14; 9 13; 10 nan; 11 15; 12 nan; 13 16; 14 17; 15 18; 16 12];
    rewardOrder{8} = [1 11; 2 12; 3 8; 4 9; 5 10; 6 14; 7 15; 8 16; 9 17; 10 nan];
    rewardOrder{11} = [1 11; 2 12; 3 8; 4 9; 5 10; 6 14; 7 15; 8 16; 9 17; 10 nan];
    pixelsPerCm =  2.7732*0.393701;

    animalNumber = str2double(projectName(4:end));
    dayNumber = str2double(dayName(4:end));

    configNumber = rewardOrder{animalNumber};
    rewardsSchema = rewardsConfig{configNumber(find(configNumber(:,1)==dayNumber),2)};
    rewardsPrevSchema = rewardsConfig{configNumber(find(configNumber(:,1)==dayNumber-1),2)};
    % doorSchema = [-5 6; -1 6; -1 10; -5 10; -5 6];


    % Get reward well positions
    animalNumber = str2double(projectName(4:end));
    dayNumber = str2double(dayName(4:end));
    schema2pixel0 = @(x) (((x-[8,8]).*[1,-1])/sqrt(58)*66+[-6 7]);
    schema2pixel = @(x) schema2pixel0(x).*[1 0.93] + [0, -5.06]; % correct for empirically observed distortions

    rewardsCurrent = schema2pixel(rewardsSchema);
    rewardsPrev = schema2pixel(rewardsPrevSchema);
    for j=1:2
        if j==1, rewards = rewardsCurrent; else, rewards = rewardsPrev; end
        xy = [behavior.position.x(:) behavior.position.y(:)];
        %         door = mean(xy(xy(:,1)<-80,:));
        door = [quantile(xy(xy(:,1)<-80,1),0.75) quantile(xy(xy(:,1)<-80,2),0.5)]; % horizontally, 75% to the right (towards the maze); vertically, the midpoint
        % order reward positions
        distance2door = sqrt(sum((mean(door,1) - rewards).^2,2));
        [~,order] = sort(distance2door);
        distance2first = sqrt(sum((mean(rewards(order==1,:),1) - rewards).^2,2));
        [~,order2] = sort(distance2first);

        % correctSequence = [mean(door,1); rewards(order==1,:); rewards(order==3,:); rewards(order==2,:); mean(door,1)];
        correctSequence = [mean(door,1); rewards(order==1,:); rewards(order2==2,:); rewards(order2==3,:); mean(door,1)];

        % allowedPoints = []; nPoints = 100;
        % for i = 1:size(correctSequence,1)-1
        %     allowedPoints = [allowedPoints; [linspace(correctSequence(i,1),correctSequence(i+1,1),nPoints)' linspace(correctSequence(i,2),correctSequence(i+1,2),nPoints)']];
        % end
        % distances = sqrt(bsxfun(@minus,xy(:,1),allowedPoints(:,1)').^2 + bsxfun(@minus,xy(:,2),allowedPoints(:,2)').^2);
        % [~,nearestNeighbour] = min(distances,[],2);
        % estimated = scores(scored(nearestNeighbour));

        % estimate distance to each trajectory:
        pt = xy; pt(:,3)= 0; d = nan(size(xy,1),size(correctSequence,1)-1);
        for i = 1:size(correctSequence,1)-1
            v1 = [correctSequence(i,:) 0]; v2 = [correctSequence(i+1,:) 0];
            q = (cross(repmat(v1 - v2,size(pt,1),1),pt - v2)) ./ norm(v1 - v2);
            d(:,i) = abs(q(:,3));
        end
        d(isnan(d)) = 0;
        d = Smooth(d,[120 0]);
        [~,trajectory] = min(d,[],2); % which of the 4 trajectories (0->1, 1->2, 2->3, 3->0) does a given point belong to
        projection = xy; l0 = nan(size(xy(:,1)));
        for i = 1:size(correctSequence,1)-1
            ok = trajectory==i;
            angle = atan2(diff(correctSequence(i:i+1,1)),diff(correctSequence(i:i+1,2)));
            correctRotated = RotateCoordinates(correctSequence(i:i+1,:),-angle);
            rotated = RotateCoordinates(xy(ok,:),-angle); rotated(:,1) = correctRotated(1,1);
            projection(ok,:) = RotateCoordinates(rotated,angle);
            rotated = (rotated(:,2) - correctRotated(1,2))./diff(correctRotated(:,2));
            rotated(rotated<0) = 0; rotated(rotated>1) = 1;
            l0(ok) = (i-1)+ rotated;
        end
        interGoalDistance = cumsum(sqrt(sum(diff(correctSequence).^2,2)));
        l = interp1([0:size(correctSequence,1)-1]',[0;interGoalDistance],l0);
        dd = sqrt(sum((projection-xy).^2,2))/pixelsPerCm;

        if j==1,
            behavior.position.currentRewards.linearized = l; behavior.position.currentRewards.distance = dd;
            behavior.position.currentRewards.projection = projection;
            behavior.position.currentRewards.boundaries = interGoalDistance;
        else
            behavior.position.previousRewards.linearized = l; behavior.position.previousRewards.distance = dd;
            behavior.position.previousRewards.projection = projection;
            behavior.position.previousRewards.boundaries = interGoalDistance;
        end
    end
    save(fullfile(basepath,[basenameFromBasepath(basepath) '.animal.behavior.mat']),'behavior');
end

MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);

try
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
catch
    SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',[]);
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
end
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

run = behavior.run;

nonrun = SubtractIntervals([0 Inf],run);

try
    try load(fullfile(basepath,'thetaCyclesTask.mat'),'cycles');
    catch load(fullfile(basepath,'thetaCyclesTrack.mat'),'cycles');
    end
catch  disp('No theta file!')
    keyboard
end
bad = InIntervals(mean(cycles,2), SubtractIntervals([0 Inf],run));
cycles(bad,:) = [];

if ~any(InIntervals(cycles,preProbeTrials)) || ~any(InIntervals(cycles,trials)),
    %     keyboard
    error
    if false
        channel = 24; channel_1_indixing = channel+1;
        nonsleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);
        intervals = nonsleep;
        for i=1:size(nonsleep)
            lfpCell{i,1} = GetAyaLFP(channel,'restrict',nonsleep(i,:));
        end
        lfp = cell2mat(lfpCell);
        [cycles,troughs] = FindThetaCycles(lfp);
        PETH(lfp,cycles(:,1),'durations',[-1 1]*2);
        PETH(lfp,cycles(:,1),'durations',[-1 1]*1); drawnow
        save('thetaCyclesTask.mat','cycles','troughs','channel','channel_1_indixing','intervals');
        clf; drawnow
        dbcont
    end
end

[spikes,regionID,regionNames,spikesCell,order] = GetAyaSpikes(pwd);

%%

% figure
clf
set(gcf,'position',[1920, 40, 1920, 960])
options = [5 10 20 25 30 40 50 Inf]; % how many cm away do we discount the positions

Qok = false(200,6); Qok(75:99,2:3) = true; Qok(101:125,4:5) = true;
Qcontrol = false(200,6); Qcontrol(75:99,4:5) = true; Qcontrol(101:125,2:3) = true;

t = behavior.timestamps(:);
try
    runLCell = {behavior.position.currentRewards.run,behavior.position.previousRewards.run};
catch
    behaviorFileChanged = true;
    for j=1:2,
        if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
            boundaries = behavior.position.currentRewards.boundaries;
        else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
            boundaries = behavior.position.previousRewards.boundaries;
        end
        ok = ~isnan(l);
        speedL = LinearVelocity([t(ok) l(ok) zeros(sum(ok),1)],5); t_0 = speedL(:,1);
        runL = t_0(FindInterval(speedL(:,2)>10));
        segments = [0; boundaries]; segments = [segments(1:end-1) segments(2:end)];
        % remove inter-segment jumps:
        [in,w] = InIntervals(t,runL);
        maxes = Accumulate(w(in),l(in),'mode','max');
        mins = Accumulate(w(in),l(in),'mode','min');
        [~,segmentMax] = InIntervals(maxes,segments);[~,segmentMin] = InIntervals(mins,segments);
        jumps = segmentMax~=segmentMin | any([segmentMax segmentMax]==0,2);
        runL(jumps,:) = [];
        nonrunL = SubtractIntervals([0 Inf],runL);
        runLCell{j} = SubtractIntervals(run,nonrunL); % running in projected space (not perpendicular to it)
    end
    behavior.position.currentRewards.run = runLCell{1};
    behavior.position.previousRewards.run = runLCell{2};
end

if behaviorFileChanged, save(fullfile(basepath,[basenameFromBasepath(basepath) '.animal.behavior.mat']),'behavior'); end

%%

% channel = 7; channel_1_indixing = channel+1;
% nonsleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);
% intervals = nonsleep;
% for i=1:size(nonsleep)
%     lfpCell{i,1} = GetAyaLFP(channel,'restrict',nonsleep(i,:));
% end
% lfp = cell2mat(lfpCell);
% [cycles,troughs] = FindThetaCycles(lfp);
% PETH(lfp,cycles(:,1),'durations',[-1 1]*2);
% clf
% PETH(lfp,cycles(:,1),'durations',[-1 1]*1); drawnow
% save('thetaCyclesTask.mat','cycles','troughs','channel','channel_1_indixing','intervals');

name = [sessionID '_OML_cheese_theta_rewards.mat'];
try
    %     error
    load(['M:\home\raly\Documents\code\new\tempFiles\' name],'scores','averages');
catch
    for kkk=optionIndex %1:length(options)
        distanceThreshold = options(kkk);
        for j=1:2
            try
                if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
                else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
                end
                closeEnough = t(FindInterval(dd<distanceThreshold));

                pos = [t,l];
                pos(:,2) = ZeroToOne(pos(:,2));
                pos1 = nani(pos);

                runL = runLCell{j};
                subplot(2,8,kkk*2+j-2);
                [splitCycles] = SplitIntervals(cycles,'nPieces',6);
                id = repmat((1:6)',length(cycles),1);
                if j==1
                    in = repelem(InIntervals(mean(cycles,2),trials) & InIntervals(mean(cycles,2),closeEnough) & InIntervals(mean(cycles,2),runL),6,1);
                else
                    in = repelem(InIntervals(mean(cycles,2),preProbeTrials) & InIntervals(mean(cycles,2),closeEnough) & InIntervals(mean(cycles,2),runL),6,1);
                end

                [estimations,actual,errors,average] = ReconstructPosition(pos1,spikes,splitCycles(in,:),'training',runL,'id',id(in));
                errorArray = reshape(errors,size(errors,1),6,[]);
                di = Shrink([nan;diff(actual)>0],6,1); di(:,2) = Shrink([nan;diff(actual)<0],6,1);
                errorArray(:,:,di(:,2)>di(:,1)) = flipud(errorArray(:,:,di(:,2)>di(:,1))); % flip cycles with more reverse movement than forward
                average = nanmean(errorArray,3);
                imagesc(repmat(average,1,2))
                try
                    if conditionType==1, title(['STIM: max distance ' num2str(distanceThreshold) 'cm'  num2str(j)]); else title(['CONTROL: max distance ' num2str(distanceThreshold) 'cm']); end
                catch
                    title(['max distance ' num2str(distanceThreshold) 'cm'])
                end

                [score,p,a,b] = FindReplayScore(average,'circular','off','nShuffles',1);
                PlotColorMap(repmat(average,1,2));
                hold on;
                plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
                plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
                set(gca,'ytick',100,'yticklabel','0','xtick','');
                PlotHVLines(100,'h','w--','linewidth',2);
                ylabel('(decoded position) - (current position)');
                % title(strrep([sessionID ' stim OFF'],'_','-'));
                if j==1
                    try
                        if conditionType==1, title(['STIM: max distance ' num2str(distanceThreshold) 'cm']); else title(['CONTROL: max distance ' num2str(distanceThreshold) 'cm']); end
                    catch
                        title(['max distance ' num2str(distanceThreshold) 'cm'])
                    end
                else
                    title(['Pre-probe < ' num2str(distanceThreshold) 'cm away'])
                end
                if kkk==1 && j==2, xlabel(basepath); end
                drawnow


                this = errorArray; this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
                this = errorArray; this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);

                scores{kkk,j} = score; averages{kkk,j} = average;
            end
        end
    end
end

scores = scores(optionIndex,:); averages = averages(optionIndex,:);

% if ~exist(fullfile('M:\home\raly\results\OML\theta',['Cheese-' sessionID '-Replay-reward-theta-sequences.fig']),'file');
% SaveFig(fullfile('M:\home\raly\results\OML\theta',['Cheese-' sessionID '-Replay-reward-theta-sequences']))
% end

% return
probeTrials = [behavior.preProbeTrials;behavior.postProbeTrials];

%% Replay

for probesOrTask = 1 % we can train the bayesian decoder on probe trials or on the task (including stim)
    for kkk=optionIndex %1:length(options)
        distanceThreshold = options(kkk);
        for j=1:2
            if probesOrTask==1
                if j==1, name = [sessionID '_OML_cheese_Cell_current_rewards_' num2str(distanceThreshold) '_cm_probeTrained.mat'];
                else name = [sessionID '_OML_cheese_Cell_prev_rewards_' num2str(distanceThreshold) '_cm_probeTrained.mat'];
                end
            elseif probesOrTask==2
                if j==1, name = [sessionID '_OML_cheese_Cell_current_rewards_' num2str(distanceThreshold) '_cm_taskTrained.mat'];
                else name = [sessionID '_OML_cheese_Cell_prev_rewards_' num2str(distanceThreshold) '_cm_taskTrained.mat'];
                end
            else
                keyboard
            end
            try
                load(['M:\home\raly\Documents\code\new\tempFiles\' name],'outputs','bursts');
                pre = InIntervals(bursts,preSleep);
                post = InIntervals(bursts,postSleep);
                r = bursts(:,[1 3]);
                durations = diff(r,[],2);
                disp([datestr(clock) ': ' name ' loaded!']);
            catch
                disp([datestr(clock) ': Computing ' name '...']);
                return
                if ~exist('bursts','var')
                    bursts = FindBursts(Restrict(spikes,sws,'shift','on'),'thresholds',[0 3],'smooth',0.01);
                    bursts(:) = Unshift(bursts(:),sws);
                    duration = diff(bursts(:,[1 3]),[],2);
                    bursts(duration>1,:) = [];
                end
                pre = InIntervals(bursts,preSleep);
                post = InIntervals(bursts,postSleep);
                r = bursts(:,[1 3]); leeway = 0.01;
                durations = diff(r,[],2);
                intervals = r; intervals0 = intervals;
                intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
                intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
                [rWindows,rID] = SplitIntervals(intervals,'pieceSize',0.02);
                if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
                else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
                end
                closeEnough = t(FindInterval(dd<distanceThreshold));

                pos = [t,l];
                pos(:,2) = ZeroToOne(pos(:,2));
                pos = nani(pos);
                runL = runLCell{j};

                if probesOrTask==1
                    [rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',Restrict(runL,probeTrials));
                elseif probesOrTask==2
                    [rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',Restrict(runL,trials));
                end
                rEstimations(isnan(rEstimations)) = 0;
                empty = (max(rEstimations)==min(rEstimations))';
                rEstimationsCell = cell(size(r,1),1);
                for i=1:size(r,1)
                    rEstimationsCell{i,1} = rEstimations(:,rID==i);
                end

                %     error('!')
                outputs = cell(size(r,1),13);
                parfor i=1:size(r,1)
                    if sum(rID==i & ~empty)>=3
                        try
                        outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
                        catch
                            keyboard
                        end
                    end
                    if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
                end
                save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos');
            end
            if j==1, outputsCurrent = outputs; else
                outputsPrev = outputs;
            end
            end
        end
end

%%
% % function [errors,g,basepath,conditionType,taskType] = BatchCanCheeseReplay(basepath,conditionType,taskType)
% 
% function [scores,averages,conditionType,taskType,basepath] = BatchCanCheeseReplay(basepath,conditionType,taskType)
% % function [errors,outputs,pre,post,durations,conditionType,taskType,basepath] = BatchCanCheeseReplay(basepath,conditionType,taskType)
% 
% %%
% options = [5 10 20 25 30 40 50 Inf];
% 
% if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
% basename = basenameFromBasepath(basepath);
% [parentFolder,dayName] = fileparts(basepath);
% [~,projectName] = fileparts(parentFolder);
% sessionID = [projectName '_' dayName];
% cd(basepath);
% 
% disp([datestr(clock) ': Starting session ' basepath '...']);
% 
% try
%     behavior = getStruct(basepath,'animal.behavior');
% catch
%     disp('No behavior file!')
%     error
% end
% try
%     trials = behavior.trials;
%     preProbeTrials = behavior.preProbeTrials;
%     if isempty(trials), error; end
% catch
%     load('cheesbTrialsDayOffset.mat')
%     trials = [cheesbTrialsDay.trials.start cheesbTrialsDay.trials.end];
%     behavior.trials = trials;
%     preProbeTrials = [cheesbTrialsDay.preProbe.start cheesbTrialsDay.preProbe.end];
%     behavior.preProbeTrials = preProbeTrials;
%     behavior.postProbeTrials = [cheesbTrialsDay.postProbe.start cheesbTrialsDay.postProbe.end];
%     save(fullfile(basepath,[basenameFromBasepath(basepath) '.animal.behavior.mat']),'behavior');
% end
% 
% % trials are sometimes wrong. Take the whole epoch
% MergePoints = getStruct(basepath,'MergePoints');
% trials = ConsolidateIntervals(MergePoints.timestamps(CountInIntervals(trials(:,1),MergePoints.timestamps)>0,:),'epsilon',1);
% % pre/post probe _epoch_ is sometimes also wrong (e.g. OML5day6). Take all running before/after trials
% preProbeTrials = [0 trials(1)];
% postProbeTrials = [trials(end) Inf];
% 
% 
% try % check if linearized positions are saved
%     l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
%     l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
%     boundaries = behavior.position.currentRewards.boundaries;
%     boundaries = behavior.position.previousRewards.boundaries;
% catch % compute linearized positions relative to current or past rewards
%     % reward configurations:
%     rewardsConfig{6} = [5,8; 9,4; 11, 12];
%     rewardsConfig{7} = [5 ,4; 11,8; 8 ,13];
%     rewardsConfig{8} = [5, 13; 8, 8; 12, 4];
%     rewardsConfig{9} = [3, 5; 10, 5; 9, 13];
%     rewardsConfig{10} = [5, 11; 13, 9; 8, 3];
%     rewardsConfig{11} = [5, 4; 8, 8; 14, 4];
%     rewardsConfig{12} = [4, 8; 12, 7; 10, 15];
%     rewardsConfig{13} = [3, 3; 8, 10; 14, 10];
%     rewardsConfig{14} = [6, 8; 11, 5; 5, 14];
%     rewardsConfig{15} = [7, 2; 13, 7; 9, 12];
%     rewardsConfig{16} = [4, 11; 13, 13; 8, 6];
%     rewardsConfig{17} = [4, 4; 12, 8; 7, 14];
%     rewardsConfig{18} = [4, 12; 8, 8; 11, 3];
%     % each OML animal is a cell. The cell contains a lookup table of [day,reward configuration]
%     rewardOrder{1} = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 11; 12 12; 13 13; 14 14; 15 15; 16 16; 17 17; 18 18; 19 nan];
%     rewardOrder{3} = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 11; 12 12; 13 13; 14 14; 15 15; 16 16; 17 17; 18 18; 19 nan];
%     rewardOrder{5} = [1 6; 2 7; 3 8; 4 9; 5 10; 6 11; 7 12; 8 13; 9 14; 10 15; 11 nan; 12 16; 13 17; 14 nan];
%     rewardOrder{7} = [1 6; 2 7; 3 8; 4 9; 5 11; 6 11; 7 12; 8 14; 9 13; 10 nan; 11 15; 12 nan; 13 16; 14 17; 15 18; 16 12];
%     rewardOrder{8} = [1 11; 2 12; 3 8; 4 9; 5 10; 6 14; 7 15; 8 16; 9 17; 10 nan];
%     rewardOrder{11} = [1 11; 2 12; 3 8; 4 9; 5 10; 6 14; 7 15; 8 16; 9 17; 10 nan];
%     pixelsPerCm =  2.7732*0.393701;
% 
%     animalNumber = str2double(projectName(4:end));
%     dayNumber = str2double(dayName(4:end));
% 
%     configNumber = rewardOrder{animalNumber};
%     rewardsSchema = rewardsConfig{configNumber(find(configNumber(:,1)==dayNumber),2)};
%     rewardsPrevSchema = rewardsConfig{configNumber(find(configNumber(:,1)==dayNumber-1),2)};
%     % doorSchema = [-5 6; -1 6; -1 10; -5 10; -5 6];
% 
% 
%     % Get reward well positions
%     animalNumber = str2double(projectName(4:end));
%     dayNumber = str2double(dayName(4:end));
%     schema2pixel0 = @(x) (((x-[8,8]).*[1,-1])/sqrt(58)*66+[-6 7]);
%     schema2pixel = @(x) schema2pixel0(x).*[1 0.93] + [0, -5.06]; % correct for empirically observed distortions
% 
%     rewardsCurrent = schema2pixel(rewardsSchema);
%     rewardsPrev = schema2pixel(rewardsPrevSchema);
%     for j=1:2
%         if j==1, rewards = rewardsCurrent; else, rewards = rewardsPrev; end
%         xy = [behavior.position.x(:) behavior.position.y(:)];
%         %         door = mean(xy(xy(:,1)<-80,:));
%         door = [quantile(xy(xy(:,1)<-80,1),0.75) quantile(xy(xy(:,1)<-80,2),0.5)]; % horizontally, 75% to the right (towards the maze); vertically, the midpoint
%         % order reward positions
%         distance2door = sqrt(sum((mean(door,1) - rewards).^2,2));
%         [~,order] = sort(distance2door);
%         distance2first = sqrt(sum((mean(rewards(order==1,:),1) - rewards).^2,2));
%         [~,order2] = sort(distance2first);
% 
%         % correctSequence = [mean(door,1); rewards(order==1,:); rewards(order==3,:); rewards(order==2,:); mean(door,1)];
%         correctSequence = [mean(door,1); rewards(order==1,:); rewards(order2==2,:); rewards(order2==3,:); mean(door,1)];
% 
%         % allowedPoints = []; nPoints = 100;
%         % for i = 1:size(correctSequence,1)-1
%         %     allowedPoints = [allowedPoints; [linspace(correctSequence(i,1),correctSequence(i+1,1),nPoints)' linspace(correctSequence(i,2),correctSequence(i+1,2),nPoints)']];
%         % end
%         % distances = sqrt(bsxfun(@minus,xy(:,1),allowedPoints(:,1)').^2 + bsxfun(@minus,xy(:,2),allowedPoints(:,2)').^2);
%         % [~,nearestNeighbour] = min(distances,[],2);
%         % estimated = scores(scored(nearestNeighbour));
% 
%         % estimate distance to each trajectory:
%         pt = xy; pt(:,3)= 0; d = nan(size(xy,1),size(correctSequence,1)-1);
%         for i = 1:size(correctSequence,1)-1
%             v1 = [correctSequence(i,:) 0]; v2 = [correctSequence(i+1,:) 0];
%             q = (cross(repmat(v1 - v2,size(pt,1),1),pt - v2)) ./ norm(v1 - v2);
%             d(:,i) = abs(q(:,3));
%         end
%         d(isnan(d)) = 0;
%         d = Smooth(d,[120 0]);
%         [~,trajectory] = min(d,[],2); % which of the 4 trajectories (0->1, 1->2, 2->3, 3->0) does a given point belong to
%         projection = xy; l0 = nan(size(xy(:,1)));
%         for i = 1:size(correctSequence,1)-1
%             ok = trajectory==i;
%             angle = atan2(diff(correctSequence(i:i+1,1)),diff(correctSequence(i:i+1,2)));
%             correctRotated = RotateCoordinates(correctSequence(i:i+1,:),-angle);
%             rotated = RotateCoordinates(xy(ok,:),-angle); rotated(:,1) = correctRotated(1,1);
%             projection(ok,:) = RotateCoordinates(rotated,angle);
%             rotated = (rotated(:,2) - correctRotated(1,2))./diff(correctRotated(:,2));
%             rotated(rotated<0) = 0; rotated(rotated>1) = 1;
%             l0(ok) = (i-1)+ rotated;
%         end
%         interGoalDistance = cumsum(sqrt(sum(diff(correctSequence).^2,2)));
%         l = interp1([0:size(correctSequence,1)-1]',[0;interGoalDistance],l0);
%         dd = sqrt(sum((projection-xy).^2,2))/pixelsPerCm;
% 
%         if j==1,
%             behavior.position.currentRewards.linearized = l; behavior.position.currentRewards.distance = dd;
%             behavior.position.currentRewards.projection = projection;
%             behavior.position.currentRewards.boundaries = interGoalDistance;
%         else
%             behavior.position.previousRewards.linearized = l; behavior.position.previousRewards.distance = dd;
%             behavior.position.previousRewards.projection = projection;
%             behavior.position.previousRewards.boundaries = interGoalDistance;
%         end
%     end
%     save(fullfile(basepath,[basenameFromBasepath(basepath) '.animal.behavior.mat']),'behavior');
% end
% 
% MergePoints = getStruct(basepath,'MergePoints');
% sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);
% 
% try
%     load(fullfile(basepath,[basename '.SleepState.states.mat']));
% catch
%     SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',[]);
%     load(fullfile(basepath,[basename '.SleepState.states.mat']));
% end
% sws = SleepState.ints.NREMstate;
% preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
% postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
% 
% run = behavior.run;
% 
% nonrun = SubtractIntervals([0 Inf],run);
% 
% try
%     try load(fullfile(basepath,'thetaCyclesTask.mat'),'cycles');
%     catch load(fullfile(basepath,'thetaCyclesTrack.mat'),'cycles');
%     end
% catch  disp('No theta file!')
%     keyboard
% end
% bad = IntervalsIntersect(cycles, SubtractIntervals([0 Inf],run));
% cycles(bad,:) = [];
% 
% if ~any(InIntervals(cycles,preProbeTrials)) || ~any(InIntervals(cycles,trials)),
%     %     keyboard
%     error
%     if false
%         channel = 17; channel_1_indixing = channel+1;
%         nonsleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);
%         intervals = nonsleep;
%         for i=1:size(nonsleep)
%             lfpCell{i,1} = GetAyaLFP(channel,'restrict',nonsleep(i,:));
%         end
%         lfp = cell2mat(lfpCell);
%         [cycles,troughs] = FindThetaCycles(lfp);
%         clf
%         PETH(lfp,cycles(:,1),'durations',[-1 1]*2);
%         PETH(lfp,cycles(:,1),'durations',[-1 1]*1); drawnow
%         save('thetaCyclesTask.mat','cycles','troughs','channel','channel_1_indixing','intervals');
%         clf; drawnow
%         dbcont
%     end
% end
% 
% [spikes,regionID,regionNames,spikesCell,order] = GetAyaSpikes(pwd);
% 
% %%
% 
% % figure
% 
% options = [5 10 20 25 30 40 50 Inf]; % how many cm away do we discount the positions
% 
% Qok = false(200,6); Qok(75:99,2:3) = true; Qok(101:125,4:5) = true;
% Qcontrol = false(200,6); Qcontrol(75:99,4:5) = true; Qcontrol(101:125,2:3) = true;
% 
% t = behavior.timestamps(:);
% for j=1:2,
%     if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
%         boundaries = behavior.position.currentRewards.boundaries;
%     else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
%         boundaries = behavior.position.previousRewards.boundaries;
%     end
%     ok = ~isnan(l);
%     speedL = LinearVelocity([t(ok) l(ok) zeros(sum(ok),1)],5); t_0 = speedL(:,1);
%     runL = t_0(FindInterval(speedL(:,2)>10));
%     segments = [0; boundaries]; segments = [segments(1:end-1) segments(2:end)];
%     % remove inter-segment jumps:
%     [in,w] = InIntervals(t,runL);
%     maxes = Accumulate(w(in),l(in),'mode','max');
%     mins = Accumulate(w(in),l(in),'mode','min');
%     [~,segmentMax] = InIntervals(maxes,segments);[~,segmentMin] = InIntervals(mins,segments);
%     jumps = segmentMax~=segmentMin | any([segmentMax segmentMax]==0,2);
%     runL(jumps,:) = [];
%     nonrunL = SubtractIntervals([0 Inf],runL);
%     runLCell{j} = SubtractIntervals(run,nonrunL); % running in projected space (not perpendicular to it)
% end
% 
% %%
% 
% name = [sessionID '_OML_cheese_theta_rewards.mat'];
% clf
% set(gcf,'position',[1920, 40, 1920, 960])
% 
% try
% %     error
%     load(['M:\home\raly\Documents\code\new\tempFiles\' name],'scores','averages');
% catch
%     for kkk=1:length(options)
%         distanceThreshold = options(kkk);
%         for j=1:2
%             try
%                 if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
%                 else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
%                 end
%                 closeEnough = t(FindInterval(dd<distanceThreshold));
% 
%                 pos = [t,l];
%                 pos(:,2) = ZeroToOne(pos(:,2));
%                 pos1 = nani(pos);
% 
%                 runL = runLCell{j};
%                 subplot(2,8,kkk*2+j-2);
%                 [splitCycles] = SplitIntervals(cycles,'nPieces',6);
%                 id = repmat((1:6)',length(cycles),1);
%                 if j==1
%                     in = repelem(InIntervals(mean(cycles,2),trials) & InIntervals(mean(cycles,2),closeEnough) & InIntervals(mean(cycles,2),runL),6,1);
%                 else
%                     in = repelem(InIntervals(mean(cycles,2),preProbeTrials) & InIntervals(mean(cycles,2),closeEnough) & InIntervals(mean(cycles,2),runL),6,1);
%                 end
% 
%                 [estimations,actual,errors,average] = ReconstructPosition(pos1,spikes,splitCycles(in,:),'training',runL,'id',id(in));
%                 errorArray = reshape(errors,size(errors,1),6,[]);
%                 di = Shrink([nan;diff(actual)>0],6,1); di(:,2) = Shrink([nan;diff(actual)<0],6,1);
%                 errorArray(:,:,di(:,2)>di(:,1)) = flipud(errorArray(:,:,di(:,2)>di(:,1))); % flip cycles with more reverse movement than forward
%                 average = nanmean(errorArray,3);
%                 imagesc(repmat(average,1,2))
%                 try
%                     if conditionType==1, title(['STIM: max distance ' num2str(distanceThreshold) 'cm'  num2str(j)]); else title(['CONTROL: max distance ' num2str(distanceThreshold) 'cm']); end
%                 catch
%                     title(['max distance ' num2str(distanceThreshold) 'cm'])
%                 end
% 
%                 [score,p,a,b] = FindReplayScore(average,'circular','off','nShuffles',1);
%                 PlotColorMap(repmat(average,1,2));
%                 hold on;
%                 plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
%                 plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
%                 set(gca,'ytick',100,'yticklabel','0','xtick','');
%                 PlotHVLines(100,'h','w--','linewidth',2);
%                 ylabel('(decoded position) - (current position)');
%                 % title(strrep([sessionID ' stim OFF'],'_','-'));
%                 if j==1
%                     try
%                         if conditionType==1, title(['STIM: max distance ' num2str(distanceThreshold) 'cm']); else title(['CONTROL: max distance ' num2str(distanceThreshold) 'cm']); end
%                     catch
%                         title(['max distance ' num2str(distanceThreshold) 'cm'])
%                     end
%                 else
%                     title(['Pre-probe < ' num2str(distanceThreshold) 'cm away'])
%                 end
%                 if kkk==1 && j==2, xlabel(basepath); end
%                 drawnow
% 
% 
%                 this = errorArray; this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
%                 this = errorArray; this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
% 
%                 scores{kkk,j} = score; averages{kkk,j} = average;
%             end
%         end
%     end
%     save(['M:\home\raly\Documents\code\new\tempFiles\' name],'scores','averages');
%     SaveFig(fullfile('M:\home\raly\results\OML\theta',['Cheese-' sessionID '-Replay-reward-theta-sequences']))
% end
% 
% % if ~exist(fullfile('M:\home\raly\results\OML\theta',['Cheese-' sessionID '-Replay-reward-theta-sequences.fig']),'file');
% % SaveFig(fullfile('M:\home\raly\results\OML\theta',['Cheese-' sessionID '-Replay-reward-theta-sequences']))
% % end
% 
% % return
% probeTrials = [behavior.preProbeTrials;behavior.postProbeTrials];
% 
% %% Replay
% for kkk=2
%     distanceThreshold = options(kkk);
%     for j=1:2
%         try
%             if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
%             else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
%             end
%             closeEnough = t(FindInterval(dd<distanceThreshold));
% 
%             pos = [t,l];
%             pos(:,2) = ZeroToOne(pos(:,2));
%             pos = nani(pos);
% 
%             runL = runLCell{j};
% 
%             % bursts = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
%             bursts = FindBursts(Restrict(spikes,sws,'shift','on'),'thresholds',[0 3],'smooth',0.01);
%             bursts(:) = Unshift(bursts(:),sws);
%             duration = diff(bursts(:,[1 3]),[],2);
%             bursts(duration>1,:) = [];
% 
%             pre = InIntervals(bursts,preSleep);
%             post = InIntervals(bursts,postSleep);
% 
%             r = bursts(:,[1 3]); leeway = 0.01;
%             durations = diff(r,[],2);
%             intervals = r; intervals0 = intervals;
%             intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
%             intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
%             [rWindows,rID] = SplitIntervals(intervals,'pieceSize',0.02);
% 
% 
%             if j==1, name = [sessionID '_OML_cheese_Cell_current_rewards_' num2str(distanceThreshold) '_cm_taskTrained.mat'];
%             else name = [sessionID '_OML_cheese_Cell_prev_rewards_' num2str(distanceThreshold) '_cm_taskTrained.mat'];
%             end
%             [rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',Restrict(runL,trials));
% 
% %             if j==1, name = [sessionID '_OML_cheese_Cell_current_rewards_' num2str(distanceThreshold) '_cm_probeTrained.mat'];
% %             else name = [sessionID '_OML_cheese_Cell_prev_rewards_' num2str(distanceThreshold) '_cm_probeTrained.mat'];
% %             end
% %             [rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',Restrict(runL,probeTrials));
%             empty = (max(rEstimations)==min(rEstimations))';
%             tic
%             rEstimationsCell = cell(size(r,1),1);
%             for i=1:size(r,1)
%                 rEstimationsCell{i,1} = rEstimations(:,rID==i);
%             end
%             try
%                 load(['M:\home\raly\Documents\code\new\tempFiles\' name],'outputs');
%                 if size(outputs,1)~=size(rEstimationsCell,1)
%                     disp(['Saved replay scores don''t fit the detected events. Recomputing...']);
%                     error(['Saved replay scores don''t fit the detected events. Recomputing...']);
%                 end
%             catch
%                 %     error('!')
%                 outputs = cell(size(r,1),13);
%                 parfor i=1:size(r,1)
%                     if sum(rID==i & ~empty)>=3
%                         outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
%                     end
%                     if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
%                 end
%                 save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos');
%             end
%             toc
%             outputsSaved = outputs; rEstimationsCellSaved = rEstimationsCell;
% 
%             %
%             durationBins = Accumulate(rID);
%             for k=1:2
%                 bad = cellfun(@isempty,outputs(:,1));
%                 if k==1
%                     ok = ~bad & pre;
%                 else
%                     ok = ~bad & post;
%                 end
%                 seshID = ones(size(ok));
%                 %         subplot(2,2,k+(condition-1)*2);
%                 scores = cell2mat(outputs(ok,1));
%                 pValues = cell2mat(outputs(ok,2));
%                 nShuffles = size(outputs{find(ok,1),5},1);
%                 shuffleID = repmat((1:nShuffles)',size(scores,1),1);
%                 shuffledScores = cell2mat(outputs(ok,5));
%                 duration = durations(ok);
%                 span = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))./durationBins(ok);
%                 slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(average,1)*4./duration; % in m/s
%                 shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(average,1)*4./repelem(duration,nShuffles);
%                 normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
%                 normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
%                 these{1,k} = [scores pValues slopes seshID(ok)];
%             end
%             %
%             clf
%             subplot(1,3,1);
%             [score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
%             PlotColorMap(Smooth(repmat(average,1,2),[2 0],'type','cc'));
%             hold on;
%             plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
%             plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
%             set(gca,'ytick',100,'yticklabel','0','xtick','');
%             PlotHVLines(100,'h','w--','linewidth',2);
%             ylabel('(decoded position) - (current position)');
%             % title(strrep([sessionID ' stim OFF'],'_','-'));
%             try
%                 if conditionType==1, title(['STIM: max distance ' num2str(distanceThreshold) 'cm']); else title(['CONTROL: max distance ' num2str(distanceThreshold) 'cm']); end
%             catch
%                 title(['max distance ' num2str(distanceThreshold) 'cm'])
%             end
%             drawnow
% 
%             %     these = saved{variety};
%             g0 = Group(these{:});
%             g = g0; % this will be normalized
%             % normalize replay
%             ok = ismember(g(:,end),[1 2]);
%             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i;
%                 normalizationFactor = [nanmedian(g0(ok,1)) diff(quantile(g0(ok,1),[0.25 0.75]))]; % subtract the median and divide by the quantile of the pre-sleep of the respective session
%                 g(ok & g(:,4)==i,1) = (g(ok & g(:,4)==i,1)-normalizationFactor(1))./normalizationFactor(2);
%             end
% 
%             subplot(1,3,2);
%             ok = g0(:,end)>0;
%             anovabar(g(ok,1),g0(ok,end),'parametric',false)
%             title(['Replay: ' strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(g(g(:,end)==1,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==2,1)))/1000)]);
%             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,n=' num2str(sum(g(:,end)==1))],['post,n=' num2str(sum(g(:,end)==2))]},'box','off');
%             ylabel('score relative to baseline');
% 
%             drawnow
% 
%             subplot(2,3,3);
%             % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
%             sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
%             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
%                 sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
%                 meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
%             end
%             isPost = g(:,end)>0;
%             anovabar(sig(isPost)-1,g(isPost,end))
%             q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
%             p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
%             text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
%             text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
%             title(['p_d=' num2str(p)]);
%             p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
%             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,p=' num2str(p(1))],['post,p=' num2str(p(2))]},'box','off');
%             set(gca,'yticklabel',get(gca,'ytick')+1);
%             set(gca,'yticklabel',get(gca,'ytick')+1);
%             ylabel('proportion replay (post/baseline)');
%             drawnow
% 
% 
%             subplot(2,3,3);
%             % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
%             sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
%             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
%                 sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
%                 meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
%             end
%             isPost = g(:,end)>0;
%             anovabar(sig(isPost)-1,g(isPost,end))
%             q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
%             p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
%             text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
%             text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
%             title(['p_d=' num2str(p)]);
%             p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
%             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,p=' num2str(p(1))],['post,p=' num2str(p(2))]},'box','off');
%             set(gca,'yticklabel',get(gca,'ytick')+1);
%             set(gca,'yticklabel',get(gca,'ytick')+1);
%             ylabel('proportion replay (post/baseline)');
%             drawnow
% 
%             subplot(2,3,6);
%             % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
%             sig = double(g(:,2)<0.05 & floor(abs(g(:,3)))>0); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
%             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
%                 sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
%                 meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
%             end
%             isPost = g(:,end)>0;
%             anovabar(sig(isPost)-1,g(isPost,end))
%             q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
%             p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
%             text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
%             text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
%             title(['p_d=' num2str(p)]);
%             p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
%             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,p=' num2str(p(1))],['post,p=' num2str(p(2))]},'box','off');
%             set(gca,'yticklabel',get(gca,'ytick')+1);
%             set(gca,'yticklabel',get(gca,'ytick')+1);
%             ylabel('proportion nonstationary (>1m/s) replay (post/baseline)');
%             drawnow
%         end
%     end
% end
% 
% % %% temp
% % for kkk=2
% %     for j=1:2
% %         try
% %             if j==1, l = behavior.position.currentRewards.linearized; dd = behavior.position.currentRewards.distance;
% %             else l = behavior.position.previousRewards.linearized; dd = behavior.position.previousRewards.distance;
% %             end
% %             closeEnough = t(FindInterval(dd<distanceThreshold));
% % 
% %             pos = [t,l];
% %             pos(:,2) = ZeroToOne(pos(:,2));
% %             pos = nani(pos);
% % 
% %             runL = runLCell{j};
% % 
% %             % bursts = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
% %             bursts = FindBursts(Restrict(spikes,sws,'shift','on'),'thresholds',[0 3],'smooth',0.01);
% %             bursts(:) = Unshift(bursts(:),sws);
% %             duration = diff(bursts(:,[1 3]),[],2);
% %             bursts(duration>1,:) = [];
% % 
% %             pre = InIntervals(bursts,preSleep);
% %             post = InIntervals(bursts,postSleep);
% % 
% %             r = bursts(:,[1 3]); leeway = 0.01;
% %             durations = diff(r,[],2);
% %             intervals = r; intervals0 = intervals;
% %             intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
% %             intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
% %             [rWindows,rID] = SplitIntervals(intervals,'pieceSize',0.02);
% % 
% %             if j==1, name = [sessionID '_OML_cheese_Cell_current_rewards_' num2str(distanceThreshold) '_cm_taskTrained.mat'];
% %             else name = [sessionID '_OML_cheese_Cell_prev_rewards_' num2str(distanceThreshold) '_cm_taskTrained.mat'];
% %             end
% %             [rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',Restrict(runL,trials));
% %             empty = (max(rEstimations)==min(rEstimations))';
% %             tic
% %             rEstimationsCell = cell(size(r,1),1);
% %             for i=1:size(r,1)
% %                 rEstimationsCell{i,1} = rEstimations(:,rID==i);
% %             end
% %             try
% %                 load(['M:\home\raly\Documents\code\new\tempFiles\' name],'outputs');
% %                 if size(outputs,1)~=size(rEstimationsCell,1)
% %                     disp(['Saved replay scores don''t fit the detected events. Recomputing...']);
% %                     error(['Saved replay scores don''t fit the detected events. Recomputing...']);
% %                 end
% %             catch
% %                 %     error('!')
% %                 outputs = cell(size(r,1),13);
% %                 parfor i=1:size(r,1)
% %                     if sum(rID==i & ~empty)>=3
% %                         outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
% %                     end
% %                     if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
% %                 end
% %                 save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos1');
% %             end
% %             toc
% %             outputsSaved = outputs; rEstimationsCellSaved = rEstimationsCell;
% % 
% %             %
% %             durationBins = Accumulate(rID);
% %             for k=1:2
% %                 bad = cellfun(@isempty,outputs(:,1));
% %                 if k==1
% %                     ok = ~bad & pre;
% %                 else
% %                     ok = ~bad & post;
% %                 end
% %                 seshID = ones(size(ok));
% %                 %         subplot(2,2,k+(condition-1)*2);
% %                 scores = cell2mat(outputs(ok,1));
% %                 pValues = cell2mat(outputs(ok,2));
% %                 nShuffles = size(outputs{find(ok,1),5},1);
% %                 shuffleID = repmat((1:nShuffles)',size(scores,1),1);
% %                 shuffledScores = cell2mat(outputs(ok,5));
% %                 duration = durations(ok);
% %                 span = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))./durationBins(ok);
% %                 slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(average,1)*4./duration; % in m/s
% %                 shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(average,1)*4./repelem(duration,nShuffles);
% %                 normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
% %                 normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
% %                 these{1,k} = [scores pValues slopes seshID(ok)];
% %             end
% %             %
% %             clf
% %             subplot(1,3,1);
% %             [score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
% %             PlotColorMap(Smooth(repmat(average,1,2),[2 0],'type','cc'));
% %             hold on;
% %             plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
% %             plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
% %             set(gca,'ytick',100,'yticklabel','0','xtick','');
% %             PlotHVLines(100,'h','w--','linewidth',2);
% %             ylabel('(decoded position) - (current position)');
% %             % title(strrep([sessionID ' stim OFF'],'_','-'));
% %             try
% %                 if conditionType==1, title(['STIM: max distance ' num2str(distanceThreshold) 'cm']); else title(['CONTROL: max distance ' num2str(distanceThreshold) 'cm']); end
% %             catch
% %                 title(['max distance ' num2str(distanceThreshold) 'cm'])
% %             end
% %             drawnow
% % 
% %             %     these = saved{variety};
% %             g0 = Group(these{:});
% %             g = g0; % this will be normalized
% %             % normalize replay
% %             ok = ismember(g(:,end),[1 2]);
% %             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i;
% %                 normalizationFactor = [nanmedian(g0(ok,1)) diff(quantile(g0(ok,1),[0.25 0.75]))]; % subtract the median and divide by the quantile of the pre-sleep of the respective session
% %                 g(ok & g(:,4)==i,1) = (g(ok & g(:,4)==i,1)-normalizationFactor(1))./normalizationFactor(2);
% %             end
% % 
% %             subplot(1,3,2);
% %             ok = g0(:,end)>0;
% %             anovabar(g(ok,1),g0(ok,end),'parametric',false)
% %             title(['Replay: ' strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(g(g(:,end)==1,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==2,1)))/1000)]);
% %             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,n=' num2str(sum(g(:,end)==1))],['post,n=' num2str(sum(g(:,end)==2))]},'box','off');
% %             ylabel('score relative to baseline');
% % 
% %             drawnow
% % 
% %             subplot(2,3,3);
% %             % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
% %             sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
% %             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
% %                 sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
% %                 meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
% %             end
% %             isPost = g(:,end)>0;
% %             anovabar(sig(isPost)-1,g(isPost,end))
% %             q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
% %             p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
% %             text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% %             text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% %             title(['p_d=' num2str(p)]);
% %             p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
% %             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,p=' num2str(p(1))],['post,p=' num2str(p(2))]},'box','off');
% %             set(gca,'yticklabel',get(gca,'ytick')+1);
% %             set(gca,'yticklabel',get(gca,'ytick')+1);
% %             ylabel('proportion replay (post/baseline)');
% %             drawnow
% % 
% % 
% %             subplot(2,3,3);
% %             % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
% %             sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
% %             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
% %                 sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
% %                 meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
% %             end
% %             isPost = g(:,end)>0;
% %             anovabar(sig(isPost)-1,g(isPost,end))
% %             q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
% %             p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
% %             text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% %             text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% %             title(['p_d=' num2str(p)]);
% %             p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
% %             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,p=' num2str(p(1))],['post,p=' num2str(p(2))]},'box','off');
% %             set(gca,'yticklabel',get(gca,'ytick')+1);
% %             set(gca,'yticklabel',get(gca,'ytick')+1);
% %             ylabel('proportion replay (post/baseline)');
% %             drawnow
% % 
% %             subplot(2,3,6);
% %             % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
% %             sig = double(g(:,2)<0.05 & floor(abs(g(:,3)))>0); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
% %             for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
% %                 sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
% %                 meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
% %             end
% %             isPost = g(:,end)>0;
% %             anovabar(sig(isPost)-1,g(isPost,end))
% %             q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
% %             p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
% %             text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% %             text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% %             title(['p_d=' num2str(p)]);
% %             p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
% %             set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['baseline,p=' num2str(p(1))],['post,p=' num2str(p(2))]},'box','off');
% %             set(gca,'yticklabel',get(gca,'ytick')+1);
% %             set(gca,'yticklabel',get(gca,'ytick')+1);
% %             ylabel('proportion nonstationary (>1m/s) replay (post/baseline)');
% %             drawnow
% %         end
% %     end
% % end
% 
% % errors = reshape(errorArray,size(errorArray,1),[]);
% % %
% % if ~exist(fullfile('M:\home\raly\results\OML\replay',['Cheese-' sessionID '-Replay-reward-barplots.fig']),'file')
% %     SaveFig(fullfile('M:\home\raly\results\OML\replay',['Cheese-' sessionID '-Replay-reward-barplots']))
% % end
% % 
% 
% 
% 
% 
% 
% 
