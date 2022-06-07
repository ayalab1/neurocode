% script_BatchCanThetaCompression

batch = StartBatch(@BatchCanThetaCompression,'OMLproject.batch');
X0 = get(batch,'UserData'); % [dPeaks,dxCorr,dTheta,cmPerTrack,basepath]

%%
figure
% X = X0; X(14:16,:) = [];
% sessionID = 'OML18';
% sessionID = 'OML19';
% sessionID = 'OLM21';
% sessionID = 'OML22';
sessionID = 'O';
X = X0;
% X = X0(14:16,:);
X = X0; X(14:16,:) = [];
 
ok = cellfun(@(x) ~isempty(strfind(x,sessionID)),X(:,end));

X(~ok,:) = []; 


[dPeaks,dxCorr,dTheta,ratioCCG,specificity] = deal(cell(1,2));
nSessions = size(X,1);
clear temp; for i=1:nSessions, temp{i,1} = reshape(ones(size(X{i,1}{1})),[],1)*i; end; seshID = cat(1,temp{:});
this = cellfun(@(x) reshape(x,[],1),cat(1,X{:,1}),'uniformoutput',0); for j=1:2,dPeaks{j} = cat(1,this{:,j}); end
this = cellfun(@(x) reshape(x,[],1),cat(1,X{:,2}),'uniformoutput',0); for j=1:2,dxCorr{j} = cat(1,this{:,j}); end
this = cellfun(@(x) reshape(x,[],1),cat(1,X{:,3}),'uniformoutput',0); for j=1:2,dTheta{j} = cat(1,this{:,j}); end
this = cellfun(@(x) reshape(x,[],1),cat(1,X{:,4}),'uniformoutput',0); for j=1:2,ratioCCG{j} = cat(1,this{:,j}); end
this = [cellfun(@(x) reshape(bsxfun(@min,x(:,1),x(:,1)'),[],1),X(:,5),'uniformoutput',0) cellfun(@(x) reshape(bsxfun(@min,x(:,2),x(:,2)'),[],1),X(:,5),'uniformoutput',0)];  for j=1:2,specificity{j} = cat(1,this{:,j}); end
% this = cellfun(@(x) reshape(x,[],1),cat(1,X{:,6}),'uniformoutput',0); for j=1:2,overlapping{j} = cat(1,this{:,j}); end
this = cellfun(@(x) reshape(x,[],1),cat(1,X{:,6}),'uniformoutput',0); for j=1:2,dCenter{j} = cat(1,this{:,j}); end

clf;
xx = [-1 1]*250; yy = [-1 1]*0.1*1000; xxx = [-1 1]*200;
nBins = [1 1]*100; smooth = 3;
[slopes,rs,ps] = deal(nan(nSessions,2)); clear maps
conditionNames = {'stim ON','stim OFF'};
for j=1:2
    subplot(3,4,j);
    data = [dPeaks{j}(:,1),dTheta{j}(:,1) ratioCCG{j}(:,1) specificity{j}(:,1) seshID];
    data(any(isnan(data),2),:) = [];
%     data = sortrows(data);
%     data(data(:,1)==0,:) =[];
%     data(data(:,4)<0.2,:) = [];
%     data = Restrict(data,xx);
    data(:,2) = data(:,2)*1000; % convert to ms
    data(data(:,3)>0.5,:) = [];
    data(:,1) = data(:,1)+randn(size(data(:,1)))/2;
    scatter(data(:,1),data(:,2),2,'filled'); colormap(jet);
%     hold on
    [r,p] = corr(data(:,1),data(:,2));
    beta = data(:,1)\data(:,2);
    title({[conditionNames{j} ': slope=' num2str(beta(1)) ', '],['r=' num2str(r) ', p=' num2str(p)]});
    xlim(xxx); ylim(yy);
%     plot(xlim,xlim*beta(1),'k--');
%     plot(xlim,xlim*(1),'k--');
    xlabel('distance from place field peaks (cm)');
    ylabel('Theta time lag (ms)');
    set(gca,'box','off','tickdir','out','fontsize',12);


    subplot(3,4,j+4);
    xBins = linspace(xx(1),xx(2),nBins(1))'; yBins = linspace(yy(1),yy(2),nBins(2))';
    [map,x,y] = DensityMap(FindClosest(xBins,data(:,1))/nBins(1),FindClosest(yBins,data(:,2))/nBins(2),'smooth',smooth,'nBins',nBins,'show','off','type','ccc'); hold on
    PlotColorMap(map/nanmean(map(:)),'x',xBins,'y',yBins);
    xlim(xxx); ylim(yy);
    %     plot(xlim,xlim*beta(1)+beta(2),'k--');
    title([conditionNames{j} ': pooled data (all pairs)']);
    xlabel('distance from place field peaks (cm)');
    ylabel('Theta time lag (ms)');
    set(gca,'box','off','tickdir','out','fontsize',12);

    subplot(3,4,j+8);
    b = Bin(data(:,1),10);
    hold all
    anovabar(data(:,2),b,'alpha',[0 0]);

    for i=1:nSessions, 
        subdata = data(data(:,end)==i,:);
        if size(subdata,1)<2, continue; end
        slopes(i,j) = subdata(:,1)\subdata(:,2);
        try [rs(i,j),ps(i,j)] = corr(subdata(:,1),subdata(:,2)); end

        try [map,x,y] = DensityMap(FindClosest(xBins,subdata(:,1))/nBins(1),FindClosest(yBins,subdata(:,2))/nBins(2),'smooth',smooth,'nBins',nBins,'show','off','type','ccc'); hold on
            maps{j}(:,:,i) = map/nanmean(map(:));
        catch
            maps{j}(1:nBins,1:nBins,i) = nan;
        end
    end
    dataCell{j} = data;
    %     subplot(2,4,j+6);
    %     hist(slopes(:,j),100);
end

slopes(ps>0.05) = nan;

subplot(3,4,7);
PlotColorMap(nanmean(maps{1},3),'x',xBins,'y',yBins)
[beta,R2,p] = CircularRegression(xBins(:),yBins(:)/1000,nanmean(maps{1},3));
title({['mean ' conditionNames{1} ' session : slope=' num2str(beta(1)*1000) ', '],['r=' num2str(r) ', p=' num2str(p)]}); 
xlim(xxx);
set(gca,'box','off','tickdir','out','fontsize',12);

subplot(3,4,8)
PlotColorMap(nanmean(maps{2},3),'x',xBins,'y',yBins)
[beta,e,p] = CircularRegression(xBins(:),yBins(:)/1000,nanmean(maps{2},3));
title({['mean ' conditionNames{2} ' session : slope=' num2str(beta(1)*1000) ', '],['r=' num2str(r) ', p=' num2str(p)]}); 
xlim(xxx);
set(gca,'box','off','tickdir','out','fontsize',12);

clims(7,8)
clims(2,5)

subplot(3,6,4); cla
g = Group(rs(:,1),rs(:,2)); g = nani(g);
anovaplot(g(:,1),g(:,2),'parametric',false);
hold on
plot(rs','.k-');
ylabel('r-coefficients per session');
set(gca,'box','off','tickdir','out','fontsize',12);
names = conditionNames; for j=1:2, names{j} = [names{j} ', p=' num2str(signrank(g(g(:,2)==j)))]; end
set(gca,'xtick',1:2,'XTickLabel',names);
title(['ranksum test: p=' num2str(ranksum(g(g(:,2)==1),g(g(:,2)==2)))]);
set(gca,'box','off','tickdir','out');

subplot(3,6,5); cla
g = Group(ps(:,1),ps(:,2));
g(isnan(g(:,1)),:) = []; g(:,1) = g(:,1)<0.05;
anovabar(g(:,1),g(:,2),'alpha',[0 0.05]);
title(['binomial test: p=' num2str(z2p(zBinomialComparison(sum(ps(:,2)<0.05),sum(~isnan(ps(:,2))),sum(ps(:,1)<0.05),sum(~isnan(ps(:,1))))))]);
ylabel('proportion of sessions with significant correlation');
names = conditionNames; for j=1:2, names{j} = [names{j} ', p=' num2str(sum(ps(:,j)<0.05)) '/' num2str(sum(~isnan(ps(:,j))))]; end
set(gca,'xtick',1:2,'XTickLabel',names);
set(gca,'box','off','tickdir','out');

subplot(3,6,6); cla
g = Group(slopes(:,1),slopes(:,2)); g = nani(g);
anovaplot(g(:,1),g(:,2),'parametric',false);
hold on
plot(slopes','.k-');
ylabel('slopes per session');
names = conditionNames; for j=1:2, names{j} = [names{j} ', p=' num2str(signrank(g(g(:,2)==j)))]; end
set(gca,'xtick',1:2,'XTickLabel',names);
title(['binomial test: p=' num2str(z2p(zBinomialComparison(sum(ps(:,2)<0.05),sum(~isnan(ps(:,2))),sum(ps(:,1)<0.05),sum(~isnan(ps(:,1))))))]);
set(gca,'box','off','tickdir','out');

drawnow
ax = axes('Units', 'normalized', 'InnerPosition', [0.01 0.5 0.08 0.4], 'Parent', gcf);
sessionNames = cellfun(@(x) x(9:end),X(:,end),'UniformOutput',0);
text(0,0.5,sessionNames);
color = get(gcf,'Color');
set(ax,'box','off','xtick',[],'ytick',[],'Color','none','XColor',color,'YColor',color,'TickDir','out')
title('Included sessions:')

% compare r-values and slopes

g = Group(dataCell{:});
labels = g(:,end);
r1 = nancorr(g(labels==1,1),g(labels==1,2));
r2 = nancorr(g(labels==2,1),g(labels==2,2));
a1 = g(labels==1,1)\g(labels==1,2);
a2 = g(labels==2,1)\g(labels==2,2);
d0 = nan(10000,2);
for i=1:length(d0), [~,order] = sort(rand(size(labels))); labels = labels(order); 
    d0(i,1) = nancorr(g(labels==2,1),g(labels==2,2)) - nancorr(g(labels==1,1),g(labels==1,2)); 
    d0(i,2) = g(labels==2,1)\g(labels==2,2) - g(labels==1,1)\g(labels==1,2); 
end

subplot(3,4,9); cla
bar([r1 r2]);
ylabel('Compression index (Pearson''s r)');
set(gca,'xtick',1:2,'XTickLabel',conditionNames);
set(gca,'box','off','tickdir','out','fontsize',10);
subplot(3,4,10); cla
hist(d0(:,1),100);
PlotHVLines(r2-r1,'v','r--','linewidth',2);
legend('shuffled','observed'); legend('location','northwest','box','off')
xlabel('difference in compression index: rOFF-rON');
title(['Permutation test: p=' num2str(Portion(d0(:,1)>=(r2-r1)))]);
set(gca,'box','off','tickdir','out','fontsize',10);

% compare slopes
drawnow
g = Group(dataCell{:});
labels = g(:,end);
r1 = g(labels==1,1)\g(labels==1,2);
r2 = g(labels==2,1)\g(labels==2,2);

subplot(3,4,11); cla
bar([a1 a2]);
ylabel('Slope (cm/ms)');
set(gca,'xtick',1:2,'XTickLabel',conditionNames);
set(gca,'box','off','tickdir','out','fontsize',10);

subplot(3,4,12); cla
hist(d0(:,2),100);
PlotHVLines(a2-a1,'v','r--','linewidth',2);
legend('shuffled','observed'); legend('location','northwest','box','off')
xlabel('difference in slope: slopeOFF-slopeON');
title(['Permutation test: p=' num2str(Portion(d0(:,2)>=(a2-a1)))]);
set(gca,'box','off','tickdir','out','fontsize',10);

% SaveFig(fullfile('M:\home\raly\results\OML\',['ThetaCompression-summary-results-excluding-678']))
