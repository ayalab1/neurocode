basepath = pwd;
kk = 0; i=1;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
for structure = 2 % 1 = hpc, 2 = pfc
    for kind = 2 % 1 for cells' firing rate, 2 for training component activity or 3, for exploration (as gauged by clicker data) component activity
        clf
        for deltaCondition = 2 % awake, all sws, Delta or No delta


            % Prepare fiigure:

            name = sessionID;
            if structure==1, name = [name ': HPC']; else, name = [name ': PFC']; end
            if kind==1, name = [name ' spikes']; else, name = [name ' 100ms']; if kind == 2, name = [name ' training']; else, name = [name ' clicker object exploration']; end; name = [name ' ICs']; end
            if deltaCondition>2, if deltaCondition==3, name = [name ' DELTA']; else, name = [name ' NO DELTA']; end; elseif deltaCondition==1,  name = [name ' awake']; else name = [name ' all SWS']; end
            figure('name',name);


            % Compute the reactivation:

            [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = BatchLoadHPCPFCData(basepath);
            if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end
            preIntervals = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
            postIntervals = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

            if deltaCondition==1, pre = SubtractIntervals(sleep(1,:), sws); post = SubtractIntervals(sleep(2:end,:), sws); end % restrict to awake periods
            try if deltaCondition>2, followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0])); if deltaCondition==3, ripples(~followed,:) = []; else, ripples(followed,:) = []; end; end
            catch, disp('no deltas/ripples'); % Issue with dividing ripples into followed and non-followed
            end

            pre = InIntervals(ripples,preIntervals); post = InIntervals(ripples,postIntervals);
            % Restrict to 1h of sws only:
            limit = Unshift(3600,preIntervals); preIntervals(preIntervals(:,1)>limit,:) = []; if preIntervals(end,2)>limit, preIntervals(end,2) = limit; end
            limit = Unshift(3600,postIntervals); postIntervals(postIntervals(:,1)>limit,:) = []; if postIntervals(end,2)>limit, postIntervals(end,2) = limit; end
            training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) contains(x,'train'),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
            try exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration); exploration = Restrict(exploration,training); end
            if structure==1, spikes = hpc; else, spikes = pfc; end
            if isempty(spikes), disp('no spikes'); continue; end
            % remove lost units
            [h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
            badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
            newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
            if any(badUnits), disp([num2str(sum(badUnits)) ' bad units found. They will be ignored.']);
                spikes(:,2) = newIDs(spikes(:,2)); spikes(spikes(:,2)==0,:) = []; % remove bad units
            end

            if kind==1 % spiking activity
                for j=1:max(spikes(:,2)) % use spiking activity
                    [q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
                    kk = kk+1;
                    qq{kk,1} = nanmean(q(pre,:),1); qq{kk,2} = nanmean(q(post,:),1); qq{kk,3} = i;
                end
            elseif kind>1 % detect components and their reactivation
                if kind==2, intervals = training; else, intervals = exploration; end
                rng(0);
                [templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
                re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
                for j=1:size(re,2)-1
                    [q,qt] = PETH(re(:,[1 1+j]),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
                    kk = kk+1;
                    qq{kk,1} = nanmean(q(pre,:),1); qq{kk,2} = nanmean(q(post,:),1); qq{kk,3} = i;
                end
            end


            % Plot the results

            clf
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
                ylim(y); PlotHVLines([0 0.2],'v','k--','linewidth',2); PlotHVLines([0],'h','k--','linewidth',2); ylim(y);
                ylabel('mean +/ s.e.m. (z-units)');
            end
            clims

            axes('position',[0.72 0.35 0.05 0.1]);
            response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
            anovaplot(response,[],'parametric',false);
            hold all
            plot(meshgrid(1:size(response,2),1:size(response,1))+randn(size(response))/50,response,'k.','markersize',15)
            set(gca,'box','off','color','none');
            drawnow

            set(gca,'TickDir','out','box','off')
            drawnow
            if strcmpi(plotType,'semplotDifference') || strcmpi(plotType,'semplot') || strcmpi(plotType,'barplots') || strcmpi(plotType,'histograms')
                EquateScales;
            end
        end
    end
end