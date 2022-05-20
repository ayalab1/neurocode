function [mPre,mPost,regionsCell,rippleRate,preferred,responses,stims,basepath] = BatchLoadV1Data(basepath)

%%
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);

session = getStruct(basepath,'session');
ripples = getStruct(basepath,'ripples');
MergePoints = getStruct(basepath,'MergePoints');
EMG = getStruct(basepath,'EMG');
immobility = EMG.timestamps(FindInterval(EMG.data<0.6)); immobility(diff(immobility,[],2)<1,:) = [];
postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
audiovisual = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'audiovisual')),MergePoints.foldernames),:),'epsilon',1);
visual = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'visual')),MergePoints.foldernames),:),'epsilon',1);
visualPre = SubtractIntervals(visual,[audiovisual(1) Inf]);
visualPost = SubtractIntervals(visual,[0 audiovisual(end)]);

[spikes,regionID,regionNames,spikesCell] = GetAyaSpikes(basepath);
% Get a single string to include in figure labels
regionListString = regionNames; regionListString(2,:) = repmat({'->'},1,size(regionListString,2)); regionListString{end} = '';
regionListString = strcat(regionListString{:});

regionNames = regionNames(:); regionID = regionID(:);
regionsCell = regionNames(regionID);

nNeurons = max(spikes(:,2));
r = ripples.timestamps;
r = Restrict(r,immobility);
pre = InIntervals(r(:,1),sleep(1,:));
post = InIntervals(r(:,1),sleep(2,:));
nBins = 501;
m = nan(nNeurons,nBins); mPre = nan(nNeurons,nBins); mPost = nan(nNeurons,nBins);
rippleRate = [sum(pre)./dur(Restrict(immobility,sleep(1,:))) sum(post)./dur(Restrict(immobility,sleep(2,:)))];

for i=1:nNeurons
    [h,ht] = PETH(spikes(spikes(:,2)==i),ripples.timestamps(:,1),'nBins',nBins,'durations',[-1 1]*5);
    m(i,:) = mean(h);
    mPre(i,:) = mean(h(pre,:));
    mPost(i,:) = mean(h(post,:));
end

%%
clf
matrix = dlmread('stimuli_synced.times');
matrix(isnan(matrix(:,1)),:) = [];
id = zeros(size(matrix(:,1)));
id(InIntervals(matrix(:,1),visualPre)) = 1;
id(InIntervals(matrix(:,1),audiovisual)) = 2;
id(InIntervals(matrix(:,1),visualPost)) = 3;
matrix(:,end) = id;

nUnits = length(spikesCell);
ok = diff([0;matrix(:,2)])~=0 & matrix(:,end)>0 & matrix(:,3)==0 & matrix(:,2)>-2;
stims = matrix(ok,:);
u = unique(stims(:,2)); %u(u==-2) = [];
for i=1:length(spikesCell)
    subplot(4,ceil(length(spikesCell)/4),i); cla
%     response = CountInIntervals(spikesCell{i},bsxfun(@plus,stims(:,1),[0 0.2]));
    response0 = CountInIntervals(spikesCell{i},bsxfun(@plus,stims(:,1),[0.065 0.165]));

    match = FindClosest(u,stims(:,2)); 
    
    cla
    response = [response0;response0;response0];
    id = [stims(:,[2 end]); stims(:,2)+180 stims(:,end); stims(:,2)+360 stims(:,end)];
    id(id==179 | id==359)=-1; id(id==-1) = -20;

    colors = {'k','b','r'};
    for k=1:3
        ok = id(:,end)==k; 
        try
        semplotBin(id(ok,1),response(ok),0,colors{k}); hold on
        catch
            keyboard
        end
        ok = stims(:,end)==k; 
        m = Accumulate(match(ok),response0(ok),'mode','mean');
        preferred(i,k) = u(findmax(m));
    end
    m = Accumulate(match,response0,'mode','mean');
    preferred(i,k+1) = u(findmax(m));
    PlotHVLines([0 180 360]+preferred(i,1),'v','k--');
    xlim([-20 360]);
    title(i);
    responses(:,i) = response0;
end


















