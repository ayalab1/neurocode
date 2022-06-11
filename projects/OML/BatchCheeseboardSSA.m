function [anglesCell,anglesCellShuffle,savedMatrix,savedMatrixShuffled,savedArray,savedArrayShuffled,m,z,basepath,conditionType,taskType] = BatchCheeseboardSSA(basepath,conditionType,taskType)

% batch = StartBatch(@BatchCheeseboardReactivation,'OMLcheese.batch');

%%
rng(0);

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
timescale = 0.015;

disp([datestr(clock) ': starting session ' basepath '...']);

cd(basepath);
MergePoints = getStruct(basepath,'MergePoints');
% ripples = getStruct(basepath,'ripples')
[spikes,regionID,regionNames,spikesCell,order] = GetAyaSpikes(basepath);
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,'cheesbTrialsDayOffset.mat'),'cheesbTrialsDay');
trials = [cheesbTrialsDay.trials.start cheesbTrialsDay.trials.end];
% postprobe = [cheesbTrialsDay.postProbe.start cheesbTrialsDay.postProbe.end];

try
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
catch
    SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',[]);
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
end

sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

% make sure the order is presleep -> task -> postsleep
preSleep = SubtractIntervals(preSleep,[trials(1) Inf]); postSleep = SubtractIntervals(postSleep,[0 trials(end)]);
assemblies = SimSpikeAssemblies(Restrict(spikes,trials,'shift','on'),timescale,sqrt(2)*erfcinv(0.05),1,0,3);
nAssemblies = size(assemblies,1);
activationIntervals = SSA_Activations(spikes,assemblies,timescale);
activations = cellfun(@(x) mean(x,2),activationIntervals,'UniformOutput',0); activations = Group(activations{:});

pre = InIntervals(activations,preSleep);
post = InIntervals(activations,postSleep);
task = InIntervals(activations,trials);

m = [Accumulate(activations(:,2),pre,'size',nAssemblies)/dur(preSleep) Accumulate(activations(:,2),task,'size',nAssemblies)/dur(trials) Accumulate(activations(:,2),post,'size',nAssemblies)/dur(postSleep)];
z = bsxfun(@rdivide,m,Accumulate(activations(:,2))/sleep(end));

%% Order stuff

triplets = SSA_SplitAssemblies(assemblies,3);
intervalsCell = {preSleep,trials,postSleep};
intervals = ConsolidateIntervals(sortrows([preSleep;trials;postSleep]));
restrictedIntervals = intervalsCell; for sub=1:length(intervalsCell), restrictedIntervals{sub}(:) = Restrict(intervalsCell{sub}(:),intervals,'shift','on'); restrictedIntervals{sub} = ConsolidateIntervals(restrictedIntervals{sub}); end
restrictedIntervals = cat(1,restrictedIntervals{:});

restrictedSpikes = Restrict(spikes,intervals,'shift','on');
restrictedSpikes(:,1) = restrictedSpikes(:,1)+randn(size(restrictedSpikes(:,1)))/1000000; % avoid ambiguities
restrictedSpikes = sortrows(restrictedSpikes);
fr = Accumulate(restrictedSpikes(:,2))/sum(diff(intervalsCell{sub},[],2));
tActivations = SSA_Activations(restrictedSpikes,triplets,timescale);

shuffled = restrictedSpikes; shuffled(:,2) = Scramble(shuffled(:,2));
tActivationsShuffled = SSA_Activations(shuffled,triplets,timescale);
for j=1:size(triplets,1)

    % Get behavior order:
    members = find(triplets(j,:));
    s = restrictedSpikes;
    clear code; code(members,1) = 1:length(members);
    code(end+1:max(s(:,2)),1) = 0;
    s(:,2) = code(s(:,2)); s(s(:,2)==0,:) = [];
    these = Restrict(tActivations{j,1},restrictedIntervals(2,:));
    sequenceList1 = GetSequenceList(s,bsxfun(@plus,these(1:2:end,:),[-1 1]*0.00001));
    sequenceList2 = GetSequenceList(s,bsxfun(@plus,these(2:2:end,:),[-1 1]*0.00001));
    sequenceList1(:,end+1:size(sequenceList2,2)) = nan; sequenceList2(:,end+1:size(sequenceList1,2)) = nan;
    sequenceList = nan(size(these,1),size(sequenceList1,2));
    sequenceList(1:2:end,:) = sequenceList1; sequenceList(2:2:end,:) = sequenceList2;
    [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,length(members),1);
    flat = reshape(matrix,[],size(matrix,3)); flat = flat(any(~isnan(flat),2),:)';
    bad = any(isnan(flat),2); %bad = bad | sum(bsxfun(@eq,sequenceList,mode(sequenceList')'),2)>1;
    flat(bad,:)=  []; sequenceList(bad,:) = [];
    [u,uu,uui] = unique(flat,'rows');
    order = sequenceList(find(uui==mode(uui),1),:); order(isnan(order)) = [];

    % Compute actual order
    s = restrictedSpikes;
    members = members(order);
%     fr_ordered(j,:) = fr(members);
    clear code;
    code(members,1) = 1:length(members);
    code(end+1:max(s(:,2)),1) = 0;
    s(:,2) = code(s(:,2)); s(s(:,2)==0,:) = [];

    [~,w] = InIntervals(mean(tActivations{j,1},2),restrictedIntervals);
    sequenceList1 = GetSequenceList(s,bsxfun(@plus,tActivations{j,1}(1:2:end,:),[-1 1]*0.00001),'keep','mid');
    sequenceList2 = GetSequenceList(s,bsxfun(@plus,tActivations{j,1}(2:2:end,:),[-1 1]*0.00001),'keep','mid');
    sequenceList1(:,end+1:size(sequenceList2,2)) = nan; sequenceList2(:,end+1:size(sequenceList1,2)) = nan;
    sequenceList = nan(size(tActivations{j,1},1),size(sequenceList1,2));
    sequenceList(1:2:end,:) = sequenceList1; sequenceList(2:2:end,:) = sequenceList2;
    [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,length(members),1);
    flat = reshape(matrix,[],size(matrix,3)); flat = flat(any(~isnan(flat),2),:)';sequenceList(bad,:) = [];
    [~,w] = InIntervals(mean(tActivations{j,1},2),restrictedIntervals);
    bad = any(isnan(flat),2);
    w(bad) = []; flat(bad,:) = [];

    % Remove nans
    sequenceList = sequenceList'; sequenceList = reshape(sequenceList(~isnan(sequenceList)),[],size(sequenceList,2))';
    [u,ui,uui] = unique(flat,'rows');
%     sequenceList(FindClosest(uui,[1:6]'),:)

    [h,ht] = Dist(1:6,[uui w],'grouped');
    %     if findmax(h(:,2))~=1,
    %         order = sequenceList(find(w==2 & uui==mode(uui(w==2)),1),:); order(isnan(order)) = [];
    %         temp{1,2} = members; temp{2,2} = order; temp{3,2} = Accumulate(uui(w==2)); temp{4,2} = uui(w==2); temp{5,2} = sequenceList(w==2,:); temp{6,2} = tActivations{j,1}(w==2,:);
    %         ['mayday ' num2str(j)]
    %         break
    %     end


    %     Dist(1:6,[uui w],'grouped'); xlim([1 6]);
    %     title(j);

    %     subplot(3,2,2);
    %     PlotColorMap(repmat(normalise(Smooth(curves(members,:),[0 2],'type','cc')),1,2));

    %     subplot(3,2,3);
    correct = 1;
    clear ps ps6 pps pps6
    %     anovabar((uui==correct)-1/6,w);
    for i=1:3, ps(i) = zBinomialComparison(sum(w==i & uui==1),sum(w==i),1/6); end
    for i=1:3, ps6(i) = zBinomialComparison(sum(w==i & uui==6),sum(w==i),1/6); end
    for i=1:3, pps(i) = sum(w==i & uui==1)/sum(w==i); end
    for i=1:3, pps6(i) = sum(w==i & uui==6)/sum(w==i); end
    %     ps(1) = -zBinomialComparison(sum(w==1 & uui==1),sum(w==1),sum(w==2 & uui==1),sum(w==2));
    %     ps(2) = -zBinomialComparison(sum(w==1 & uui==1),sum(w==1),sum(w==3 & uui==1),sum(w==3));
    %     ps(3) = zBinomialComparison(sum(w==3 & uui==1),sum(w==3),sum(w==2 & uui==1),sum(w==2));
    %     title(ps);

    %     subplot(3,2,4);
    sequenceList1 = GetSequenceList(s,bsxfun(@plus,tActivations{j,1}(1:2:end,:),[-1 1]*0.001),'repetitions','on');
    sequenceList2 = GetSequenceList(s,bsxfun(@plus,tActivations{j,1}(2:2:end,:),[-1 1]*0.001),'repetitions','on');
    sequenceList1(:,end+1:size(sequenceList2,2)) = nan; sequenceList2(:,end+1:size(sequenceList1,2)) = nan;
    sequenceList = nan(size(tActivations{j,1},1),size(sequenceList1,2));
    sequenceList(1:2:end,:) = sequenceList1; sequenceList(2:2:end,:) = sequenceList2;
    [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,length(members),1);
%     [~,ref_ibeforej,ref_count] = GetBiasMatrix(1:length(members),length(members),1);
    ref_ibeforej = nanmean(ibeforej(:,:,w==2),3); ref_count = nanmean(count(:,:,w==2),3);
    c = BiasCosine(ibeforej,count,ref_ibeforej,ref_count); c(bad) = [];

    angles = acos(c);
    % find a second reference for left/right
    ref_ibeforej = ibeforej(:,:,FindClosest(c,0)); ref_count = count(:,:,FindClosest(c,0));
    c2 = BiasCosine(ibeforej,count,ref_ibeforej,ref_count); c(bad) = [];
    angles(c2<0) = -angles(c2<0);
    anglesCell{j,1} = [angles w];

    [h2,ht] = Dist(min([numel(unique(c)) 100]),[c w],'grouped');

    %     plot(ht,h2);
    %     title(ht*h2);
    %     xlim([-1 1]);

    saved(j,:) = ht*h2;
    [h2,ht] = Dist(linspace(-1,1,4),[c w],'grouped');
    saved_z(j,:) = ps;
    saved_z6(j,:) = ps6;
    saved_p(j,:) = pps; saved_p6(j,:) = pps6;
    saved_h1(:,:,j) = h;
    saved_h2(:,:,j) = h2;
    c = corr(h);
    saved_c(j,:) = [c(1,2) c(1,3) c(2,3)];

    % Compare to shuffled order
    members = find(triplets(j,:));
    s = shuffled;
    clear code; code(members,1) = 1:length(members);
    code(end+1:max(s(:,2)),1) = 0;
    s(:,2) = code(s(:,2)); s(s(:,2)==0,:) = [];
    these = Restrict(tActivationsShuffled{j,1},restrictedIntervals(2,:));
    sequenceList1 = GetSequenceList(s,bsxfun(@plus,these(1:2:end,:),[-1 1]*0.00001));
    sequenceList2 = GetSequenceList(s,bsxfun(@plus,these(2:2:end,:),[-1 1]*0.00001));
    sequenceList1(:,end+1:size(sequenceList2,2)) = nan; sequenceList2(:,end+1:size(sequenceList1,2)) = nan;
    sequenceList = nan(size(these,1),size(sequenceList1,2));
    sequenceList(1:2:end,:) = sequenceList1; sequenceList(2:2:end,:) = sequenceList2;
    [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,length(members),1);
    flat = reshape(matrix,[],size(matrix,3)); flat = flat(any(~isnan(flat),2),:)';
    bad = any(isnan(flat),2); %bad = bad | sum(bsxfun(@eq,sequenceList,mode(sequenceList')'),2)>1;
    flat(bad,:)=  []; sequenceList(bad,:) = [];
    [u,uu,uui] = unique(flat,'rows');
    order = sequenceList(find(uui==mode(uui),1),:); order(isnan(order)) = [];

    s = shuffled;
    members = members(order);
%     fr_ordered(j,:) = fr(members);
    clear code;
    code(members,1) = 1:length(members);
    code(end+1:max(s(:,2)),1) = 0;
    s(:,2) = code(s(:,2)); s(s(:,2)==0,:) = [];
    [~,w] = InIntervals(mean(tActivationsShuffled{j,1},2),restrictedIntervals);
    sequenceList1 = GetSequenceList(s,bsxfun(@plus,tActivationsShuffled{j,1}(1:2:end,:),[-1 1]*0.00001),'keep','mid');
    sequenceList2 = GetSequenceList(s,bsxfun(@plus,tActivationsShuffled{j,1}(2:2:end,:),[-1 1]*0.00001),'keep','mid');
    sequenceList1(:,end+1:size(sequenceList2,2)) = nan; sequenceList2(:,end+1:size(sequenceList1,2)) = nan;
    sequenceList = nan(size(tActivationsShuffled{j,1},1),size(sequenceList1,2));
    sequenceList(1:2:end,:) = sequenceList1; sequenceList(2:2:end,:) = sequenceList2;
    [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,length(members),1);
    flat = reshape(matrix,[],size(matrix,3)); flat = flat(any(~isnan(flat),2),:)';sequenceList(bad,:) = [];
    [~,w] = InIntervals(mean(tActivationsShuffled{j,1},2),restrictedIntervals);
    bad = any(isnan(flat),2);
    w(bad) = []; flat(bad,:) = [];
    sequenceList = sequenceList'; sequenceList = reshape(sequenceList(~isnan(sequenceList)),[],size(sequenceList,2))';
    [u,ui,uui] = unique(flat,'rows');
    [h,ht] = Dist(1:6,[uui w],'grouped');
    clear ps ps6 pps pps6
    for i=1:3, ps(i) = zBinomialComparison(sum(w==i & uui==1),sum(w==i),1/6); end
    for i=1:3, ps6(i) = zBinomialComparison(sum(w==i & uui==6),sum(w==i),1/6); end
    for i=1:3, pps(i) = sum(w==i & uui==1)/sum(w==i); end
    for i=1:3, pps6(i) = sum(w==i & uui==6)/sum(w==i); end
    sequenceList1 = GetSequenceList(s,bsxfun(@plus,tActivationsShuffled{j,1}(1:2:end,:),[-1 1]*0.001),'repetitions','on');
    sequenceList2 = GetSequenceList(s,bsxfun(@plus,tActivationsShuffled{j,1}(2:2:end,:),[-1 1]*0.001),'repetitions','on');
    sequenceList1(:,end+1:size(sequenceList2,2)) = nan; sequenceList2(:,end+1:size(sequenceList1,2)) = nan;
    sequenceList = nan(size(tActivationsShuffled{j,1},1),size(sequenceList1,2));
    sequenceList(1:2:end,:) = sequenceList1; sequenceList(2:2:end,:) = sequenceList2;
    [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,length(members),1);
    [~,ref_ibeforej,ref_count] = GetBiasMatrix(1:length(members),length(members),1);
    c = BiasCosine(ibeforej,count,ref_ibeforej,ref_count); c(bad) = [];
     angles = acos(c);
    % find a second reference for left/right
    ref_ibeforej = ibeforej(:,:,FindClosest(c,0)); ref_count = count(:,:,FindClosest(c,0));
    c2 = BiasCosine(ibeforej,count,ref_ibeforej,ref_count); c(bad) = [];
    angles(c2<0) = -angles(c2<0);
    anglesCellShuffle{j,1} = [angles w];
    [h2,ht] = Dist(min([numel(unique(c)) 100]),[c w],'grouped');
    saved0(j,:) = ht*h2;
    [h2,ht] = Dist(linspace(-1,1,4),[c w],'grouped');
    saved_z0(j,:) = ps;
    saved_z60(j,:) = ps6;
    saved_p0(j,:) = pps; 
    saved_p60(j,:) = pps6;
    saved_h10(:,:,j) = h;
    saved_h20(:,:,j) = h2;
    c = corr(h);
    saved_c0(j,:) = [c(1,2) c(1,3) c(2,3)];    
end

savedMatrix = [saved saved_p saved_p6 saved_c saved_z saved_z6];
savedArray = cat(1,saved_h1,saved_h2);
savedMatrixShuffled = [saved0 saved_p0 saved_p60 saved_c0 saved_z0 saved_z60];
savedArrayShuffled = cat(1,saved_h10,saved_h20);



















