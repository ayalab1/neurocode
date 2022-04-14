% batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
% batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
% batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
% Xo = get(batchO,'UserData'); % opto
% Xs = get(batchS,'UserData'); % short
% Xl = get(batchL,'UserData'); % long
%
% XX = {Xs,Xl,Xo};
% write = false;
%%
% figure

plotTypes = {'anovabars','semplot','barplots','postColorMap','multiple','histograms','semplotDifference'};
deltaConditionNames = {'awake ripples','all sws ripples','ripples followed by deltas','sws ripples not followed by deltas'}; conditionNames = {'short','extended','opto'};
for iii = 1
    %     clf
    plotType = 'semplot';
    %     plotType = plotTypes{iii};
    colors = {'k','r','m'};
    colorsConditions = {[5 145 157]/255,[255 175 51]/255,[128 0 128]/255};
    for structure = 2 % 1 = hpc, 2 = pfc
        for kind = 2 % 1 for cells' firing rate, 2 for training component activity or 3, for exploration (as gauged by clicker data) component activity
            clf
            for deltaCondition = 2:4% awake, all sws, Delta or No delta
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
                        training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
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
                                %                                                     qqs{kk,1} = q;
                            end
                        elseif kind>1 % detect components and their reactivation
                            if kind==2, intervals = training; else, intervals = exploration; end
                            rng(0);
                            [templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
                            re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
                            for j=1:size(re,2)-1
                                [q,qt] = PETH(re(:,[1 1+j]),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
                                kk = kk+1;
                                qq{kk,1} = nanmean(q(pre,:),1);
                                qq{kk,2} = nanmean(q(post,:),1);
                                qq{kk,3} = i;
                                %                             qqs{kk,1} = q;
                            end
                        end
                        if i==nSessions, qq(kk+1:end,:) = []; end
                    end

                    if condition==1, name = 'Short training'; elseif condition==2, name = 'Extended training'; else name = 'Opto prolonged ripples'; end
                    if structure==1
                        name = [name ': HPC'];
                    else
                        name = [name ': PFC'];
                    end
                    if kind==1
                        name = [name ' spikes'];
                    else
                        name = [name ' 100ms'];
                        if kind == 2
                            name = [name ' training'];
                        else
                            name = [name ' clicker object exploration'];
                        end
                        name = [name ' ICs'];
                    end
                    if deltaCondition>2
                        if deltaCondition==3
                            name = [name ' DELTA'];
                        else
                            name = [name ' NO DELTA'];
                        end
                    elseif deltaCondition==1,  name = [name ' awake'];
                    else name = [name ' all SWS'];
                    end

                    %
                    nComponents = Accumulate(cell2mat(qq(:,end)));
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

                    subplot(3,4,(condition-1)*4+deltaCondition);

                    if strcmpi(plotType,'barplots')
                        matrices = {z1,z2,z2-z1};
                        response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
                        kruskalbar(response);
                        if structure==1, ylabel('Median response 0-100ms (z-units)');
                        else, ylabel('Median response 0-200ms (z-units)');
                        end
                        title({name,['p=' num2str(signrank(response(:,2)-response(:,1)))]});
                        set(gca,'xticklabel',{'pre','post'});
                    end

                    if strcmpi(plotType,'histograms') || strcmpi(plotType,'anovabars')
                        if structure==1, response = [mean(matrices{1}(:,InIntervals(qt,[0 0.1])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.1])),2)];
                        else, response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
                        end
                        d = diff(response,[],2);
                        dd{deltaCondition,condition} = d;
                        if strcmpi(plotType,'histograms')
                            subplot(3,4,(condition-1)*4+deltaCondition); cla
                            hist(d,linspace(-5,5,50));
                            handle = get(gca,'Children');
                            set(handle,'FaceColor',colorsConditions{condition},'EdgeAlpha',0);
                            PlotHVLines(0,'v','k--','linewidth',2);
                            if structure==1, ylabel('post-pre response 0-100ms (z-units)');
                            else, ylabel('post-pre response 0-200ms (z-units)');
                            end
                            title({name, ['p=' num2str(signrank(d))]});
                            xlabel(['difference, median =' num2str(nanmedian(d))]);
                        end
                    end

                    if strcmpi(plotType,'semplot')
                        j=1; handle = semplot(qt,matrices{j},colors{j},smooth(end));
                        if j==1,
                            handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
                        end
                        ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
                        if j==1,legend([handle handle2],{'pre','post'}); else legend(handle,'difference'); end
                        xlabel('time from ripple start (s)');
                        title(name);
                    end

                    if strcmpi(plotType,'semplotDifference')
                        handle = semplot(qt,matrices{2}-matrices{1},colorsConditions{condition},smooth(end)); y=ylim;
                        ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
                        legend(handle,'post-pre');
                        xlabel('time from ripple start (s)');
                        title(name);
                    end


                    if strcmpi(plotType,'postcolormap')
                        j=2;
                        PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
                        set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
                        PlotHVLines(lines,'h','w--','linewidth',2);
                        title([name ', post'])
                        PlotHVLines(0,'v','k--','linewidth',2);
                    end

                    if strcmpi(plotType,'multiple')
                        clf
                        set(gcf,'position',[1 1 1920 1024])
                        matrices = {z1,z2,z2-z1};
                        titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
                        colors = {'k','r','m'};
                        for j=1:3
                            subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
                            PlotHVLines(lines,'h','w--','linewidth',2);
                            set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
                            usedClims(j,:) = clim;

                            title(titles{j});
                            subplot(2,3,3+j);
                            if j==2
                                response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
                                kruskalbar(response);
                                ylabel('Median response 0-500ms (zz-units)');
                                title(['p=' num2str(signrank(response(:,2)-response(:,1)))]);
                                set(gca,'xticklabel',{'pre','post'});
                            else
                                handle = semplot(qt,matrices{j},colors{j},smooth(end)); y=ylim;
                                if j==1, handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
                                    legend([handle handle2],{'pre','post'});
                                end
                                ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
                                %                             if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
                                %                                 if h>0,
                                %                                     intervals = (FindInterval(stats.below));
                                %                                     if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
                                %                                         plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
                                %                                     end
                                %                                     intervals = (FindInterval(stats.above));
                                %                                     if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
                                %                                         plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
                                %                                     end
                                %                                 end
                                %                             end
                                ylabel('mean +/ s.e.m. (z-units)');
                            end
                        end
                        clims([max(usedClims(:,1)) min(usedClims(:,2))])
                        if write
                            SaveFig(fullfile('M:\home\raly\results\PFC\reactivation',strrep(strrep(name,':',''),' ','_')))
                        end
                    end
                    set(gca,'TickDir','out','box','off')
                    drawnow
                    if strcmpi(plotType,'semplotDifference') || strcmpi(plotType,'semplot') || strcmpi(plotType,'barplots') || strcmpi(plotType,'histograms')
                        EquateScales;
                    end
                end
            end
        end

        if write
            if structure==1
                name = ['HPC'];
            else
                name = ['PFC'];
            end
            if kind==1
                name = [name ' spikes'];
            else
                name = [name ' 100ms'];
                if kind == 2
                    name = [name ' training'];
                else
                    name = [name ' clicker object exploration'];
                end
                name = [name ' ICs'];
            end

            if strcmpi(plotType,'barplots')
                SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_delta_no_delta_barplots']));
            end

            if strcmpi(plotType,'histograms')
                SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_delta_no_delta_differences_histograms']));
                g = Group(dd{:}); gg = g; gg(:,2) = ceil(g(:,2)/4); gg(:,3) = rem(g(:,2)-1,4)+1;
                for i=1:4,
                    subplot(2,2,i); cla; hold off;
                    ok = gg(:,end)==i;
                    means = Accumulate(gg(ok,2),g(ok,1))./Accumulate(gg(ok,2),~isnan(g(ok,1)));
                    anovabar(g(ok,1),gg(ok,2)); hold on;
                    handle1 = bar(1,means(1)); handle2 = bar(2,means(2)); handle3 = bar(3,means(3));
                    set(handle1,'facecolor',colorsConditions{1}); set(handle2,'facecolor',colorsConditions{2}); set(handle3,'facecolor',colorsConditions{3});
                    title(deltaConditionNames{i});
                    ylabel({['mean ' name],['difference (z-units)']});
                    set(gca,'xticklabel',conditionNames,'TickDir','out','box','off','fontsize',12);
                end
                EquateScales
                drawnow
                SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_delta_no_delta_differences_anovabars']));
            end

            if strcmpi(plotType,'semplot')
                SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_semplot_delta_no_delta']));
            end

            if strcmpi(plotType,'semplotDifference')
                SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_semplotDifferences_delta_no_delta']));
            end

            if strcmpi(plotType,'postcolormap')
                clims
                SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_post_colormaps_delta_no_delta']));
            end
        end
    end
end


%%


plotTypes = {'anovabars','semplot','barplots','postColorMap','multiple','histograms','semplotDifference'};
deltaConditionNames = {'awake ripples','all sws ripples','ripples followed by deltas','sws ripples not followed by deltas'}; conditionNames = {'short','extended','opto'};
for iii = 1
    %     clf
    plotType = 'semplot';
    %     plotType = plotTypes{iii};
    colors = {'k','r','m'};
    colorsConditions = {[5 145 157]/255,[255 175 51]/255,[128 0 128]/255};
    for structure = 2 % 1 = hpc, 2 = pfc
        for kind = 2 % 1 for cells' firing rate, 2 for training component activity or 3, for exploration (as gauged by clicker data) component activity
            clf
            for condition = 1:3  % for each of the three conditions (short training, long training, opto prolongation)
                X = XX{condition};
                condition
                qq = {}; sessionIDs = {}; kk = 0; nSessions = size(X,1);
                for i=1:size(X,1)
                    [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
                    if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end;
                    if isempty(deltas), continue; end
                    [parentFolder,basename] = fileparts(basepath);
                    [~,projectName] = fileparts(parentFolder);
                    sessionID = [projectName '_' basename];
                    sessionIDs{i,1} = sessionID;
                    preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
                    postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

                    ripples = Restrict(ripples,sws);
                    deltas = Restrict(deltas,sws);

                    pre = InIntervals(deltas,preSleep);
                    post = InIntervals(deltas,postSleep);
                    % Restrict to 1h of sws only:
                    limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
                    limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
                    training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
                    exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
                    exploration = Restrict(exploration,training);
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
                            [q,qt] = PETH(spikes(spikes(:,2)==j),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
                            kk = kk+1;
                            qq{kk,1} = nanmean(q(pre,:),1);
                            qq{kk,2} = nanmean(q(post,:),1);
                            qq{kk,3} = i;
                            %                                                     qqs{kk,1} = q;
                        end
                    elseif kind>1 % detect components and their reactivation
                        if kind==2, intervals = training; else, intervals = exploration; end
                        rng(0);
                        [templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
                        re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
                        for j=1:size(re,2)-1
                            [q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
                            kk = kk+1;
                            qq{kk,1} = nanmean(q(pre,:),1);
                            qq{kk,2} = nanmean(q(post,:),1);
                            qq{kk,3} = i;
                            %                             qqs{kk,1} = q;
                        end
                    end
                    if i==nSessions, qq(kk+1:end,:) = []; end
                end

                if condition==1, name = 'Short training'; elseif condition==2, name = 'Extended training'; else name = 'Opto prolonged ripples'; end
                if structure==1
                    name = [name ': HPC'];
                else
                    name = [name ': PFC'];
                end
                if kind==1
                    name = [name ' spikes'];
                else
                    name = [name ' 100ms'];
                    if kind == 2
                        name = [name ' training'];
                    else
                        name = [name ' clicker object exploration'];
                    end
                    name = [name ' ICs'];
                end

                %
                nComponents = Accumulate(cell2mat(qq(:,end)));
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

                subplot(3,4,(condition-1)*4+deltaCondition);

%                 if strcmpi(plotType,'barplots')
%                     matrices = {z1,z2,z2-z1};
%                     response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
%                     kruskalbar(response);
%                     if structure==1, ylabel('Median response 0-100ms (z-units)');
%                     else, ylabel('Median response 0-200ms (z-units)');
%                     end
%                     title({name,['p=' num2str(signrank(response(:,2)-response(:,1)))]});
%                     set(gca,'xticklabel',{'pre','post'});
%                 end
% 
%                 if strcmpi(plotType,'histograms') || strcmpi(plotType,'anovabars')
%                     if structure==1, response = [mean(matrices{1}(:,InIntervals(qt,[0 0.1])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.1])),2)];
%                     else, response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
%                     end
%                     d = diff(response,[],2);
%                     dd{deltaCondition,condition} = d;
%                     if strcmpi(plotType,'histograms')
%                         subplot(3,4,(condition-1)*4+deltaCondition); cla
%                         hist(d,linspace(-5,5,50));
%                         handle = get(gca,'Children');
%                         set(handle,'FaceColor',colorsConditions{condition},'EdgeAlpha',0);
%                         PlotHVLines(0,'v','k--','linewidth',2);
%                         if structure==1, ylabel('post-pre response 0-100ms (z-units)');
%                         else, ylabel('post-pre response 0-200ms (z-units)');
%                         end
%                         title({name, ['p=' num2str(signrank(d))]});
%                         xlabel(['difference, median =' num2str(nanmedian(d))]);
%                     end
%                 end

                if strcmpi(plotType,'semplot')
                    subplot(3,4,(0)*4 + condition);
                    j=1; handle = semplot(qt,matrices{j},colors{j},smooth(end));
                    if j==1,
                        handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
                    end
                    ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
                    if j==1,legend([handle handle2],{'pre','post'}); else legend(handle,'difference'); end
                    xlabel('time from delta peak (s)');
                    title(name);

                    subplot(3,4,(1)*4 + condition);
                    handle = semplot(qt,matrices{2}-matrices{1},colorsConditions{condition},smooth(end)); y=ylim;
                    ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
                    xlabel('time from delta peak (s)');
                    legend(handle,'post-pre');
                    title(name);
                end


                if strcmpi(plotType,'postcolormap')
                    j=2;
                    PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
                    set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
                    PlotHVLines(lines,'h','w--','linewidth',2);
                    title([name ', post'])
                    PlotHVLines(0,'v','k--','linewidth',2);
                end

                if strcmpi(plotType,'multiple')
                    clf
                    set(gcf,'position',[1 1 1920 1024])
                    matrices = {z1,z2,z2-z1};
                    titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
                    colors = {'k','r','m'};
                    for j=1:3
                        subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
                        PlotHVLines(lines,'h','w--','linewidth',2);
                        set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
                        usedClims(j,:) = clim;

                        title(titles{j});
                        subplot(2,3,3+j);
                        if j==2
                            response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
                            kruskalbar(response);
                            ylabel('Median response 0-500ms (zz-units)');
                            title(['p=' num2str(signrank(response(:,2)-response(:,1)))]);
                            set(gca,'xticklabel',{'pre','post'});
                        else
                            handle = semplot(qt,matrices{j},colors{j},smooth(end)); y=ylim;
                            if j==1, handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
                                legend([handle handle2],{'pre','post'});
                            end
                            ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
                            %                             if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
                            %                                 if h>0,
                            %                                     intervals = (FindInterval(stats.below));
                            %                                     if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
                            %                                         plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
                            %                                     end
                            %                                     intervals = (FindInterval(stats.above));
                            %                                     if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
                            %                                         plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
                            %                                     end
                            %                                 end
                            %                             end
                            ylabel('mean +/ s.e.m. (z-units)');
                        end
                    end
                    clims([max(usedClims(:,1)) min(usedClims(:,2))])
                    if write
                        SaveFig(fullfile('M:\home\raly\results\PFC\reactivation',strrep(strrep(name,':',''),' ','_')))
                    end
                end
                set(gca,'TickDir','out','box','off')
                drawnow
                if strcmpi(plotType,'semplotDifference') || strcmpi(plotType,'semplot') || strcmpi(plotType,'barplots') || strcmpi(plotType,'histograms')
                    EquateScales;
                end
            end
        end
    end
end

if write
    if structure==1
        name = ['HPC'];
    else
        name = ['PFC'];
    end
    if kind==1
        name = [name ' spikes'];
    else
        name = [name ' 100ms'];
        if kind == 2
            name = [name ' training'];
        else
            name = [name ' clicker object exploration'];
        end
        name = [name ' ICs'];
    end

    if strcmpi(plotType,'barplots')
        SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_delta_no_delta_barplots']));
    end

    if strcmpi(plotType,'histograms')
        SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_delta_no_delta_differences_histograms']));
        g = Group(dd{:}); gg = g; gg(:,2) = ceil(g(:,2)/4); gg(:,3) = rem(g(:,2)-1,4)+1;
        for i=1:4,
            subplot(2,2,i); cla; hold off;
            ok = gg(:,end)==i;
            means = Accumulate(gg(ok,2),g(ok,1))./Accumulate(gg(ok,2),~isnan(g(ok,1)));
            anovabar(g(ok,1),gg(ok,2)); hold on;
            handle1 = bar(1,means(1)); handle2 = bar(2,means(2)); handle3 = bar(3,means(3));
            set(handle1,'facecolor',colorsConditions{1}); set(handle2,'facecolor',colorsConditions{2}); set(handle3,'facecolor',colorsConditions{3});
            title(deltaConditionNames{i});
            ylabel({['mean ' name],['difference (z-units)']});
            set(gca,'xticklabel',conditionNames,'TickDir','out','box','off','fontsize',12);
        end
        EquateScales
        drawnow
        drawnow
        SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_delta_centered_differences_anovabars']));
    end

    if strcmpi(plotType,'semplot')
        SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_semplot_delta_centered']));
    end

    if strcmpi(plotType,'semplotDifference')
        SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_semplotDifferences_delta_centered']));
    end

    if strcmpi(plotType,'postcolormap')
        clims
        SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\' strrep(strrep(name,':',''),' ','_') '_post_colormaps_delta_centered']));
    end
end
