batch = StartBatch(@BatchLoadV1Data,'AudioVisual.batch');
% X = get(batch,'UserData');

%% 
X = get(batch,'UserData');
% X = X(2:end,:);
% X = X(1,:);

regions = cat(1,X{:,3});
sesh = repelem((1:size(X,1))',cellfun(@(x) size(x,1),X(:,3)));

regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
    regionCell(strcmp(regions,regionNames{i})) = i;
end

mPre = cell2mat(X(:,1));
mPost = cell2mat(X(:,2));
mDiff = mPost-mPre;
bad = all(isnan(mDiff),2);
sesh(bad) = []; regionCell(bad) = []; mPre(bad,:) = []; mPost(bad,:) = []; mDiff(bad,:) = [];

z = [mPre mPost];
zPre = (mPre - nanmean(z,2))./nanstd(z,[],2);
zPost = (mPost - nanmean(z,2))./nanstd(z,[],2);
zDiff = zPost - zPre;

zzPre = nanzscore(mPre,[],2);
zzPost = nanzscore(mPost,[],2);
% zzPre = zBaseline(mPre,ht>-2);
% zzPost = zBaseline(mPost,ht>-2);
zzDiff = zzPost-zzPre;
ht = linspace(-5,5,size(mPre,2));

figure(3)
clf
% datas = {mPre,mPost,mDiff,zPre,zPost,zDiff};
datas = {zPre,zPost,zDiff,zzPre,zzPost,zzDiff};
for i=1:6,
    data = datas{i};
    subplot(2,3,i);
    %     if i==2 | i==5, subplot(2,3,i-1); end
    titl = '';
    %         semplot(ht,data(regionCell==1,:),'r'); titl = [titl regionNames{1} ', '];
    %     semplot(ht,data(regionCell==2,:),[1 0.5 0]); titl = [titl regionNames{2} ', '];
    %     semplot(ht,data(regionCell==3,:),'g'); titl = [titl regionNames{3} ', '];
    %     semplot(ht,data(regionCell==4,:),'b'); titl = [titl regionNames{4} ', '];
    semplot(ht,data(regionCell==1,:),'k'); titl = [titl regionNames{1} ', '];
    PlotHVLines(0,'v','k--'); PlotHVLines(0,'h','k--');
    xlim([-1 1]*2);
    title(titl);
end
EquateScales

%%
clf
datas = {mPre,mPost,mDiff,zPre,zPost,zDiff};
datas = {zPre,zPost,zDiff,zzPre,zzPost,zzDiff};
for i=1:6,
    data = datas{i};
    subplot(2,3,i);
    titl = '';
    presented = []; ns = [];
    for k=1
        presented = [presented; data(regionCell==k,:)]; titl = [titl regionNames{k} ', ']; ns = [ns; sum(regionCell==k)];
    end
    PlotColorMap(Smooth(presented,[0 2]),'x',ht);
    PlotHVLines(0,'v','w--','linewidth',2);
    PlotHVLines(cumsum(ns)+0.5,'h','w','linewidth',2);
    title(titl);
    if i==1, c = clim; end
    if ismember(i,[1 2 4 5]), clim(c); end
    if rem(i,3)==0, clim([-1 1]*3); end
    xlim([-2 2]);
end

%%
clf
for i=1:3,
    data = datas{i+3};
    subplot(2,3,i);
    titl = '';
    presented = []; ns = [];
    for k=2
        presented = [presented; data(regionCell==k,:)]; titl = [titl regionNames{k} ', ']; ns = [ns; sum(regionCell==k)];
    end
    PlotColorMap(Smooth(presented,[0 2]),'x',ht);
    PlotHVLines(0,'v','w--','linewidth',2);
%     PlotHVLines(cumsum(ns)+0.5,'h','w','linewidth',2);
%     title(titl);
    if i==1, c = clim; end
    if ismember(i,[1 2 4 5]), clim(c); end
    if rem(i,3)==0, clim([-1 1]*3); end
    xlim([-2 2]);
    ylabel('CA1 units');
    set(gca,'fontsize',12,'tickdir','out','box','off');
    subplot(2,3,i+3);
    %     if i==2 | i==5, subplot(2,3,i-1); end
    titl = '';
    %         semplot(ht,data(regionCell==1,:),'r'); titl = [titl regionNames{1} ', '];
    %     semplot(ht,data(regionCell==2,:),[1 0.5 0]); titl = [titl regionNames{2} ', '];
    %     semplot(ht,data(regionCell==3,:),'g'); titl = [titl regionNames{3} ', '];
    %     semplot(ht,data(regionCell==4,:),'b'); titl = [titl regionNames{4} ', '];
    semplot(ht,data(regionCell==k,:),'k'); titl = [titl regionNames{1} ', '];
    PlotHVLines(0,'v','k--'); PlotHVLines(0,'h','k--');
    xlim([-1 1]*2);
    ylabel('mean +/- sem firing rate (z-units)');
    xlabel('time from ripple start (s)');
    set(gca,'fontsize',12,'tickdir','out','box','off');
%     title(titl);
end
% SaveFig('M:\home\raly\results\V1\Pooled\V1_later_days_ripple_response_zscored_separately');
% SaveFig('M:\home\raly\results\V1\Pooled\V1_day1_CA2_ripple_response_zscored_separately');

%%

datas = {zPre,zPost,zDiff,zzPre,zzPost,zzDiff};
for i=1:2,
    data = datas{(i)*3};
    subplot(1,2,i);
    titl = '';
    presented = []; ns = [];
    for k=1:5
        presented = [presented; data(regionCell==k,:)]; titl = [titl regionNames{k} ', ']; ns = [ns; sum(regionCell==k)];
    end
    PlotColorMap(presented,'x',ht);
    PlotHVLines(0,'v','w--','linewidth',2);
    PlotHVLines(cumsum(ns)+0.5,'h','w','linewidth',2);
    title(titl);
    clim([-1 1]*min(abs(clim)));
    xlim([-1 1]*0.5);
end

%% V1 analysis by audiovisual training day
figure;
data = zzDiff;
colors = flipud(jet(5));
k = 1; % for V1
clf
for j=1:6,
    data = datas{j}(regionCell==k,:);
    subplot(2,3,j);
    ns = Accumulate(sesh(regionCell==k));
    PlotColorMap(Smooth(data,[0 2]),'x',ht);
    PlotHVLines(0,'v','w--','linewidth',2);
    PlotHVLines(cumsum(ns)+0.5,'h','w','linewidth',2);
    xlim([-1 1]*2);
end

for j=1:6
    subplot(2,3,j);
    if j==1, c = clim*1.5; end
    if ismember(j,[1 2 4 5]), clim(c); end
    if rem(j,3)==0, clim([-1 1]*min(abs(clim))); ColorMap(gca,[0 0 1],[1 1 1],[1 0 0]); end
end


%% Code testing if the curve has shifted

for k=1:size(X,1)
    [responses,stims] = X{k,6:7};

    % clf
    % smooth = [1 0 1]*0;
    % for i=1:size(responses,2)
    %     subplot(4,ceil(size(responses,2)/4),i); cla
    %     %     response = CountInIntervals(spikesCell{i},bsxfun(@plus,stims(:,1),[0 0.2]));
    %     response0 = responses(:,i);
    %
    %     match = FindClosest(u,stims(:,2));
    %
    %     cla
    %     response = [response0;response0;response0];
    %     id = [stims(:,[2 end]); stims(:,2)+180 stims(:,end); stims(:,2)+360 stims(:,end)];
    %     id(id==179 | id==359)=-1; id(id==-1) = -20;
    %
    %     colors = {'k','b','r'};
    %     for k=1:3
    %         ok = id(:,end)==k;
    %         try
    %             semplotBin(id(ok,1),response(ok),0,colors{k},smooth(k)); hold on
    %         catch
    %             keyboard
    %         end
    %         ok = stims(:,end)==k;
    %         m = Accumulate(match(ok),response0(ok),'mode','mean');
    %         preferred(i,k) = u(findmax(m));
    %     end
    %     m = Accumulate(match,response0,'mode','mean');
    %     preferred(i,k+1) = u(findmax(m));
    %     PlotHVLines([0 180 360]+preferred(i,1),'v','k--');
    %     xlim([-20 360]);
    %     title(i);
    % end

    b = -1 + 2*ismember(stims(:,2),[0 90]) + ismember(stims(:,2),[45 135]);
    [m1,m2,preference] = deal(nan(size(response,2),2));
    for i=1:size(responses,2),
        ok1 = b==1 & stims(:,end)==1;
        ok0 = b==0 & stims(:,end)==1;
        m1(i,1) = nanmean(responses(ok1,i));
        m2(i,1) = nanmean(responses(ok0,i));
        [~,p] = ttest2(responses(ok1,i),responses(ok0,i));
        preference(i,1) = sign(nanmean(responses(ok1,i))-nanmean(responses(ok0,i))) * p2z(p);

        ok1 = b==1 & stims(:,end)==3;
        ok0 = b==0 & stims(:,end)==3;
        m1(i,2) = nanmean(responses(ok1,i));
        m2(i,2) = nanmean(responses(ok0,i));
        [~,p] = ttest2(responses(ok1,i),responses(ok0,i));
        preference(i,2) = sign(nanmean(responses(ok1,i))-nanmean(responses(ok0,i))) * p2z(p);
    end

    d = (m1-m2)./(m1+m2);

    D{k,1} = d;
    D{k,2} = preference;
end

%%
clf
set(gcf,'position',[600 400 1000 500]);
preference = cell2mat(D(:,2));
subplot(1,2,1); 
anovabar(preference,[],'parametric',false); 
set(gca,'xtick',1:2,'xticklabel',{'before','after audiovisual'});
ylabel('preference for 0 or 90 degrees vs 45 or 135 degrees');
axis square
set(gca,'fontsize',12,'tickdir','out','box','off');
title(['p=' num2str(signrank(diff(preference,[],2),0,'tail','left'))]);
subplot(1,2,2);
b = bootstrp(100000,@(x) mean(diff(x,[],2)),preference);
[h,ht] = hist(b,100); 
ci = quantile(b,[0.025 0.975]);
colors = get(gca,'ColorOrder');
handle = fill(Smooth(h,2),ht,colors(1,:));
set(handle,'facealpha',0.5,'EdgeColor','none');
hold on
plot([0 0],ci,'k','linewidth',3);
plot(0,median(b),'k.','markersize',30);
PlotHVLines(0,'h','k--');
set(gca,'box','off','xtick',[]);
axis square
set(gca,'fontsize',12,'tickdir','out','box','off');
set(gca,'xtick',[]);
ylabel('difference (preference z-units)');
title({'Bootstrapped increase (post - pre)','of preference for 0 and 90 degree stimuli'});
xlabel(['p=' num2str(mean(q<0))]);











