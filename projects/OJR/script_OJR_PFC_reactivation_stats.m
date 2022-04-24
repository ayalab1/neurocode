batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
Xo = get(batchO,'UserData'); % opto
Xs = get(batchS,'UserData'); % short
Xl = get(batchL,'UserData'); % long

XX = {Xs,Xl,Xo};
write = false;
%%
% figure

deltaConditionNames = {'awake ripples','all sws ripples','ripples followed by deltas','sws ripples not followed by deltas'}; conditionNames = {'short','extended','opto'};

colors = {'k','r','m'};
colorsConditions = {[5 145 157]/255,[255 175 51]/255,[128 0 128]/255};
clear saved savedQ
for structure = 1 % 1 = hpc, 2 = pfc
    for kind = 1 % 1 for cells' firing rate, 2 for training component activity or 3, for exploration (as gauged by clicker data) component activity
        clf
        for deltaCondition = 2%:4% awake, all sws, Delta or No delta
            for condition = 1:3  % for each of the three conditions (short training, long training, opto prolongation)
                X = XX{condition};
                condition
                qq = {}; sessionIDs = {}; kk = 0; nSessions = size(X,1);
                for i=1:size(X,1)
                    [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
                    if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end;
                    [parentFolder,basename] = fileparts(basepath);
                    [~,projectName] = fileparts(parentFolder);
                    sessionID = [projectName '_' basename];
                    sessionIDs{i,1} = sessionID;
                    preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
                    postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

                    if deltaCondition==1
                        preSleep = SubtractIntervals(sleep(1,:), sws);
                        postSleep = SubtractIntervals(sleep(2:end,:), sws);
                    end
                    try

                        if deltaCondition>2
                            followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
                            if deltaCondition==3
                                ripples(~followed,:) = [];
                            else
                                ripples(followed,:) = [];
                            end
                        end
                    catch
                        continue;
                    end

                    pre = InIntervals(ripples,preSleep);
                    post = InIntervals(ripples,postSleep);
                    % Restrict to 1h of sws only:
                    limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
                    limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
                    training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) contains(x,'train'),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
                    try exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
                        exploration = Restrict(exploration,training); end
                    if structure==1, spikes = hpc; else, spikes = pfc; end
                    if isempty(spikes), continue; end
                    % remove lost units
                    [h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
                    badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
                    newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
                    spikes(:,2) = newIDs(spikes(:,2));
                    spikes(spikes(:,2)==0,:) = []; % remove bad units

                    if kind==1
                        for j=1:max(spikes(:,2)) % use spiking activity
                            [q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
                            kk = kk+1;
                            qq{kk,1} = nanmean(q(pre,:),1);
                            qq{kk,2} = nanmean(q(post,:),1);
                            qq{kk,3} = i;
                            qq{kk,4} = sum(spikes(:,2)==j)/MergePoints.timestamps(end);
                            %                                                     qqs{kk,1} = q;
                        end
                    elseif kind>1 % detect components and their reactivation
                        if kind==2, intervals = training; else, intervals = exploration; end
                        rng(0);
                        [templates,correlations,weights,variance] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1,'tracyWidom',true,'mode','ica');
                        re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
                        for j=1:size(re,2)-1
                            [q,qt] = PETH(re(:,[1 1+j]),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
                            kk = kk+1;
                            qq{kk,1} = nanmean(q(pre,:),1);
                            qq{kk,2} = nanmean(q(post,:),1);
                            qq{kk,3} = i;
                            qq{kk,4} = variance(j);
                            %                             qqs{kk,1} = q;
                        end
                    end
                    if i==nSessions, qq(kk+1:end,:) = []; end
                end

                %
                nComponents = Accumulate(cell2mat(qq(:,3)));
                addSlash = @(x) strrep(x,'_','\_');
                sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
                lines = cumsum(nComponents(nComponents>0))+0.5;
                q = (cell2mat(qq(:,1:2)));
                %         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
                z = nanzscore(q,[],2);
                z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
                smooth = [0 2];

                z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
                matrices = {z1,z2,z2-z1};


                if condition==1, name = 'Short training'; elseif condition==2, name = 'Extended training'; else name = 'Opto prolonged ripples'; end
                if structure==1, name = [name ': HPC']; else, name = [name ': PFC']; end
                if kind==1, name = [name ' spikes']; else, name = [name ' 100ms']; if kind == 2, name = [name ' training']; else, name = [name ' clicker object exploration']; end; name = [name ' ICs']; end
                if deltaCondition>2, if deltaCondition==3, name = [name ' DELTA']; else, name = [name ' NO DELTA']; end; elseif deltaCondition==1,  name = [name ' awake']; else name = [name ' all SWS']; end
                figure('name',name);
                set(gcf,'position',[1 1 1920 1024])

                matrices = {z1,z2,z2-z1};
                titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
                colors = {'k','r','m'};
                for j=1:3
                    subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
                    PlotHVLines(lines,'h','w--','linewidth',2);
                    set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);

                    title(titles{j});
                    subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
                    if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
                    ylim(y); PlotHVLines(intervalOfInterest,'v','k--','linewidth',2); PlotHVLines([0],'h','k--','linewidth',2); ylim(y);
                    ylabel('mean +/ s.e.m. (z-units)');
                end
                clims

                axes('position',[0.72 0.35 0.05 0.1]);
                response = [mean(matrices{1}(:,InIntervals(qt,intervalOfInterest)),2) mean(matrices{2}(:,InIntervals(qt,intervalOfInterest)),2)];
                anovaplot(response,[],'parametric',false);
                set(gca,'xtick',1:2,'xticklabel',{'pre','post'});
                hold all
                plot(meshgrid(1:size(response,2),1:size(response,1))+randn(size(response))/50,response,'k.','markersize',15)
                set(gca,'box','off','color','none');
                drawnow

                set(gca,'TickDir','out','box','off')
                saved(condition,:) = matrices;
                savedQ(condition,:) = {cell2mat(qq(:,1)),cell2mat(qq(:,2)),cell2mat(qq(:,3)),cell2mat(qq(:,4))};
                drawnow

                SaveFig(fullfile('M:\home\raly\results\PFC\reactivation',strrep(strrep(['Summary ' name ' around ripples'],' ','_'),':','')));
            end
        end
    end
end


%%
figure

for i=1:3, savedQ{i,5} = i*ones(size(savedQ{i,4})); end
v = cell2mat(savedQ(:,4:5));

z = nanzscore(cell2mat(savedQ(:,1:2)),[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
zd = z2-z1;
zz1 = zscore(z1,[],2); zz2 = zscore(z2,[],2);
zzd = zz2-zz1;
intervalOfInterest = [0 0.1];

clf

% Data normalized together
subplot(2,5,1);
PlotColorMap(sortby(z1,v(:,2)*100+v(:,1)),'x',qt);
PlotHVLines(CumSum(Accumulate(v(:,2)))+0.5,'h','w--','linewidth',2);
PlotHVLines(intervalOfInterest,'v','k--','linewidth',1); 
ylabel('zscore pre/post PETH together');
subplot(2,5,2);
PlotColorMap(sortby(z2,v(:,2)*100+v(:,1)),'x',qt);
PlotHVLines(CumSum(Accumulate(v(:,2)))+0.5,'h','w--','linewidth',2);
PlotHVLines(intervalOfInterest,'v','k--','linewidth',2); 
clims(1,2)
subplot(2,5,4);
PlotColorMap(sortby(zd,v(:,2)*100+v(:,1)),'x',qt);
PlotHVLines(CumSum(Accumulate(v(:,2)))+0.5,'h','w--','linewidth',2);
PlotHVLines(intervalOfInterest,'v','k--','linewidth',1); 
clim([-1 1]*mean(abs(clim)));
subplot(2,5,5);
for i=1:3, ok = v(:,2)==i; semplot(qt,sortby([zd(ok,:)],v(ok,2)*1000+v(ok,1)),colorsConditions{i}); drawnow; pause(0.1); end
PlotHVLines(intervalOfInterest,'v','k--','linewidth',2); PlotHVLines([0],'h','k--','linewidth',2);
subplot(2,5,3);
m1 = cell2mat(savedQ(:,1)); m2 = cell2mat(savedQ(:,2));
response = [nanmean(m1(:,InIntervals(qt,intervalOfInterest)),2) nanmean(m2(:,InIntervals(qt,intervalOfInterest)),2)];
d = @(x) (x(:,2)-x(:,1));
anovaplot(d(response),v(:,2),'parametric',false);
PlotHVLines([0],'h','k--','linewidth',1); 
hold all
plot(v(:,2)+randn(size(v(:,2)))/20,d(response),'k.','markersize',20);
set(gca,'xtick',1:3,'xticklabel',{'short','long','opto'});
ylabel('[0 100]ms response difference');



% each PETH zscored individually
subplot(2,5,6);
PlotColorMap(sortby(zz1,v(:,2)*100+v(:,1)),'x',qt);
PlotHVLines(CumSum(Accumulate(v(:,2)))+0.5,'h','w--','linewidth',2);
PlotHVLines(intervalOfInterest,'v','k--','linewidth',1); 
ylabel('zscore pre/post PETH individually');
subplot(2,5,7);
PlotColorMap(sortby(zz2,v(:,2)*100+v(:,1)),'x',qt);
PlotHVLines(CumSum(Accumulate(v(:,2)))+0.5,'h','w--','linewidth',2);
PlotHVLines(intervalOfInterest,'v','k--','linewidth',1); 
clims(3,4)
subplot(2,5,9);
PlotColorMap(sortby(zzd,v(:,2)*100+v(:,1)),'x',qt);
PlotHVLines(CumSum(Accumulate(v(:,2)))+0.5,'h','w--','linewidth',2);
PlotHVLines(intervalOfInterest,'v','k--','linewidth',1); 
clim([-1 1]*mean(abs(clim)));

subplot(2,5,10);
for i=1:3, ok = v(:,2)==i; semplot(qt,sortby([zzd(ok,:)],v(ok,2)*1000+v(ok,1)),colorsConditions{i}); end
PlotHVLines(intervalOfInterest,'v','k--','linewidth',2); PlotHVLines([0],'h','k--','linewidth',2);

for i=[1:2 4:7 9:10],subplot(2,5,i); xlim([-1 1]*0.5); end

subplot(2,5,8);
response = [nanmean(zz1(:,InIntervals(qt,intervalOfInterest)),2) nanmean(zz2(:,InIntervals(qt,intervalOfInterest)),2)];
d = @(x) (x(:,2)-x(:,1));
anovaplot(d(response),v(:,2),'parametric',false);
PlotHVLines([0],'h','k--','linewidth',1); 
hold all
plot(v(:,2)+randn(size(v(:,2)))/20,d(response),'k.','markersize',20);
set(gca,'xtick',1:3,'xticklabel',{'short','long','opto'});
ylabel('[0 100]ms response difference');

% SaveFig(fullfile('M:\home\raly\results\PFC\reactivation',strrep(strrep('Summary PFC 50ms training ICs all SWS around ripples all conditions',' ','_'),':','')));
















