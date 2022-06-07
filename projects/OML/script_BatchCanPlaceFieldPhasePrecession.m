batch = StartBatch(@BatchCanPlaceFieldPhasePrecession,'OMLproject.batch');
X0 = get(batch,'UserData');

%%
% figure;
X = X0;
for i=1:size(X,1), for j=1:size(X{i,1},1), X{i,1}{j,1}(:,end+1) = i; end; end
for i=1:size(X,1), for j=1:size(X{i,1},1), X{i,1}{j,2}(:,end+1) = i; end; end

% sessionID = 'OML18';
% sessionID = 'OML19';
% sessionID = 'OLM21';
% sessionID = 'OML22';
sessionID = 'O';
X = X0;
% X = X0(14:16,:);
X = X0; X(14:16,:) = [];
% 
ok = cellfun(@(x) ~isempty(strfind(x,sessionID)),X(:,end));

X(~ok,:) = []; 

clf
points = cat(1,X{:,1});
names = {'stim ON','stim OFF'};
% these = cell2mat(points(:,1));  hmhm = Accumulate(these(:,end),these(:,1)==0,'mode','mean');
% these = cell2mat(points(:,2));  hmhm(:,2) = Accumulate(these(:,end),these(:,1)==0,'mode','mean');

smoothed = cat(1,X{:,2});
nBins = size(smoothed{find(~cellfun(@isempty,smoothed(:,1)),1),1});
% nBins = [50 20];
smooth = [1 0.5].*[50 20];
x = linspace(0,1,nBins(1));
y = linspace(0,4*pi,nBins(2)*2);

for j=1:2,
    for plotKind = 1:3
        subplot(2,3,(j-1)*3+plotKind);
        these = cell2mat(points(:,j)); these = these(:,1:2);
%         these(these(:,1)<0.01,:) = []; these(these(:,1)>0.99,:) = [];
        doubled = [these; these(:,1) these(:,2)+2*pi];
        switch plotKind
            case 1
                [beta,c,p] = CircularRegression(these(:,1),these(:,2),'slope',0);
                PlotXY(doubled(1:5:end,:),'.','markersize',1);
                title([names{j} ': pooled data'],['slope = ' num2str(beta(1)) ', intercept' num2str(beta(2)) ', r = ' num2str(c) ', p = ' num2str(p)]);
                ylim([0 4*pi]);
            case 2
                DensityMap(doubled(:,1),doubled(:,2),'smooth',smooth,'nBins',nBins,'type','clc');
                clabel('Number of points per bin');
                title({[names{j} ': pooled phase precession density plot']});
            case 3
                q = cat(3,smoothed{:,j});
                for i=1:size(q,3),
                    this = q(:,:,i);
                    this = Smooth(this,smooth,'type','lc');
                    this(:) = zscore(this(:));
                    q(:,:,i) = this;
                end
                subplot(2,3,(j-1)*3+3);
                PlotColorMap(Smooth(repmat(nanmean(q,3),2,1),0,'type','cc'),'x',x,'y',y);
%                 clim([0.5 1.3]);
                clim([-1 1]*0.3)
                clabel('Normalized density');
                title({[names{j} ': average normalized'],'phase precession density plots for all place fields'});

%                 saved{j,1} = Smooth(repmat(nanmean(q,3),2,1),smooth,'type','cc');
                saved{j,1} = q;
        end
        PlotHVLines(pi:pi:3*pi,'h','w:','linewidth',1);
        if plotKind<3
            plot(xlim,xlim*beta(1)+beta(2)+0*pi,'k--','linewidth',2);
            plot(xlim,xlim*beta(1)+beta(2)+2*pi,'k--','linewidth',2);
            plot(xlim,xlim*beta(1)+beta(2)+4*pi,'k--','linewidth',2);
        end
        xlabel('Location within place field');
        ylabel('theta phase (rad.)');
        yticks = round(linspace(0,4*pi,5)*100)/100;
        set(gca,'ytick',yticks,'YTickLabel',repmat({'    '},size(yticks)),'xtick',0:0.5:1);
        text(-0.05*ones(size(yticks)),yticks+0.4,{'0','\pi','2\pi','3\pi','4\pi'},'horizontal','center','vertical','top','fontsize',12);
        set(gca,'box','off','tickdir','out','fontsize',12);
        drawnow
    end
end

drawnow
ax = axes('Units', 'normalized', 'InnerPosition', [0.02 0.5 0.08 0.4], 'Parent', gcf);
sessionNames = cellfun(@(x) x(9:end),X(:,end),'UniformOutput',0);
text(0,0.5,sessionNames);
color = get(gcf,'Color');
set(ax,'box','off','xtick',[],'ytick',[],'Color','none','XColor',color,'YColor',color,'TickDir','out')
title('Included sessions:')

stats = cat(1,X{:,3});
g = Group(cell2mat(stats(:,1)),cell2mat(stats(:,2)));
ax = axes('Units', 'normalized', 'InnerPosition', [0.02 0.11 0.08 0.3], 'Parent', gcf);
counts = [Accumulate(g(:,end),g(:,3)<0.05 & g(:,4)<0)';Accumulate(g(:,end),g(:,3)>-1)']';
ylim([0 100]); hold on;
CountBar(counts,1);
title(['p=' num2str(z2p(zBinomialComparison(counts(2,1),counts(2,2),counts(1),counts(1,2))))]);
ylabel('Percent place fields with significant phase precession');
for j=1:2, theseNames{j} = [names{j} ': ' num2str(counts(j,1)) '/' num2str(counts(j,2))]; end
set(gca,'xtick',1:2,'xticklabel',theseNames);

% SaveFig(fullfile('M:\home\raly\results\OML\',['Phase-precession-summary-results']))
