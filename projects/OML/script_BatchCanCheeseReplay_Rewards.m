batch = StartBatch(@BatchCanCheeseReplay_Rewards,'OMLcheese.batch');

X0 = get(batch,'UserData');
%%
conditionX = 0;
X = X0(cell2mat(X0(:,9))==conditionX,:);
X(cellfun(@sum,X(:,5))<100 | cellfun(@sum,X(:,6))<100,:) = []; % at least 100 events

figure(conditionX+1+2);
sessionID = 'O';
ok = cellfun(@(x) ~isempty(strfind(x,sessionID)),X0(:,end));

X(~ok,:) = [];  %X(14:16,:) = [];
e = cat(1,X{:,2});
eon = cat(2,e{:,1});
eoff = cat(2,e{:,2});

ns = [cellfun(@(x) size(x{1},1), X(:,1)) cellfun(@(x) size(x{2},1), X(:,1))];
mon = squeeze(nansum(reshape(bsxfun(@times,eon,repelem(ns(:,1),6)'),200,6,[]),3))/sum(ns(:,1));
moff = squeeze(nansum(reshape(bsxfun(@times,eoff,repelem(ns(:,2),6)'),200,6,[]),3))/sum(ns(:,2));

scores = cat(1,X{:,1});
% figure('name','18');
clf
subplot(2,6,1);
average = mon;
[score,p,a,b] = FindReplayScore(average,'circular','off','nShuffles',1);
PlotColorMap(repmat(average,1,2));
hold on; plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep(['Theta: current correct path'],'_','-'));

subplot(2,6,2);
average = moff;
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
PlotColorMap(repmat(average,1,2));
hold on;
plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep(['Theta: previous correct path'],'_','-'));
clims([0.0045 0.008]);
drawnow

subplot(2,6,7); cla
sessionNames = cellfun(@(x) x(9:end),X(:,end),'UniformOutput',0);
text(0,0.5,sessionNames);
set(gca,'xtick',[],'ytick',[])
title('sessions included');

subplot(2,6,8);
onScore = cat(1,scores{:,1});
offScore = cat(1,scores{:,2});
g0 = Group(onScore,offScore); g0(isnan(g0(:,2)),:) = [];
g =[-diff(g0(:,1:2)./sum(g0(:,1:2),2),[],2),g0(:,3)];
anovaplot(-diff(g0(:,1:2)./sum(g0(:,1:2),2),[],2),g0(:,3),'parametric','off');
p = [signrank(g(g(:,2)==1)) signrank(g(g(:,2)==2))];
set(gca,'xtick',1:2,'xticklabel',{['ON, p=' num2str(p(1))],['OFF, p=' num2str(p(2))]});
p = ranksum(g(g(:,2)==1),g(g(:,2)==2));
title(['Theta: quadrant scores: ranksum p=' num2str(p)]);

% Replay
outputsOn = cat(1,X{:,3}); outputsOff = cat(1,X{:,4});
pre = cat(1,X{:,5}); post = cat(1,X{:,6}); durations = cat(1,X{:,7});
seshID = repelem((1:size(X,1))',cellfun(@length,X(:,6)),1);

for variety = 1
    switch variety
        case 1 % pool as bursts
            for condition = 1:2
                if condition==1
                    outputs = outputsOn;
                    bad = cellfun(@isempty,outputs(:,1));
                    bad = bad;% | ~(reOn);
                else
                    outputs = outputsOff;
                    bad = cellfun(@isempty,outputs(:,1));
                    bad = bad;% | ~(reOff);
                end
                for k=1:2
                    if k==1
                        ok = ~bad & pre;
                    else
                        ok = ~bad & post;
                    end
                    %         subplot(2,2,k+(condition-1)*2);
                    scores = cell2mat(outputs(ok,1));
                    pValues = cell2mat(outputs(ok,2));
                    nShuffles = size(outputs{find(ok,1),5},1);
                    shuffleID = repmat((1:nShuffles)',size(scores,1),1);
                    shuffledScores = cell2mat(outputs(ok,5));
                    duration = durations(ok);
                    slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(average,1)*4./duration; % in m/s
                    shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(average,1)*4./repelem(duration,nShuffles);
                    normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
                    normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
                    these{condition,k} = [scores pValues slopes seshID(ok)];
                end
            end
        case 2 % keet better fit only:
            reOn = false(size(outputsOn,1),1);
            reOff = false(size(outputsOn,1),1);
            for i=1:size(outputsOn,1), try reOn(i,1) = outputsOn{i,1}>outputsOff{i,1}; end; end
            for i=1:size(outputsOn,1), try reOff(i,1) = outputsOn{i,1}<outputsOff{i,1}; end; end
            for condition = 1:2
                if condition==1
                    outputs = outputsOn;
                    bad = cellfun(@isempty,outputs(:,1));
                    bad = bad | ~(reOn);
                else
                    outputs = outputsOff;
                    bad = cellfun(@isempty,outputs(:,1));
                    bad = bad | ~(reOff);
                end
                for k=1:2
                    if k==1
                        ok = ~bad & pre;
                    else
                        ok = ~bad & post;
                    end
                    %         subplot(2,2,k+(condition-1)*2);
                    scores = cell2mat(outputs(ok,1));
                    pValues = cell2mat(outputs(ok,2));
                    nShuffles = size(outputs{find(ok,1),5},1);
                    shuffleID = repmat((1:nShuffles)',size(scores,1),1);
                    shuffledScores = cell2mat(outputs(ok,5));
                    duration = durations(ok);
                    slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(average,1)*4./duration; % in m/s
                    shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(average,1)*4./repelem(duration,nShuffles);
                    normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
                    normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
                    these{condition,k} = [scores pValues slopes seshID(ok)];
                end
            end
    end
    saved{variety} = these;
end

% %%
% figure

for variety = 1
    these = saved{variety};
    g0 = Group(these{:});
    g = g0; % this will be normalized
    % normalize ON replay
    on = ismember(g(:,end),[1 3]);
    for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i;
        normalizationFactor = [nanmedian(g0(ok,1)) diff(quantile(g0(ok,1),[0.25 0.75]))]; % subtract the median and divide by the quantile of the pre-sleep of the respective session
        g(on & g(:,4)==i,1) = (g(on & g(:,4)==i,1)-normalizationFactor(1))./normalizationFactor(2);
    end
    % normalize OFF replay
    off = ismember(g(:,end),[2 4]);
    for i=1:max(seshID), ok = g(:,end)==2 & g(:,4)==i;
        normalizationFactor = [nanmedian(g0(ok,1)) diff(quantile(g0(ok,1),[0.25 0.75]))]; % subtract the median and divide by the quantile of the pre-sleep of the respective session
        g(off & g(:,4)==i,1) = (g(off & g(:,4)==i,1)-normalizationFactor(1))./normalizationFactor(2);
    end

    subplot(1,3,(variety-1)*3+2);
    ok = g0(:,end)>2;
    anovabar(g(ok,1),g0(ok,end)-2,'parametric',false)
    if variety==1
        title(['Replay: ' strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
    elseif variety==2
        title(['Replay: better fit (current vs prev): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
    else
        title(['Replay: decoded together (sorted posthoc): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
    end
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'current','previous'},'box','off');
    ylabel('score relative to baseline');
    ylim([0 0.4]);

    drawnow

    subplot(1,6,(variety-1)*6+5);
    % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
    sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; on = ismember(g(:,end),[1 3]);   off = ismember(g(:,end),[2 4]);
    for i=1:max(seshID),
        ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok));
        preOn(i,variety) = expected; postOn(i,variety) = nanmean(sig(g(:,end)==3 & g(:,4)==i));
        sig(g(:,4)==i & on) = sig(g(:,4)==i & on)./expected;
        ok = g(:,end)==2 & g(:,4)==i; expected = nanmean(sig(ok));
        preOff(i,variety) = expected; postOff(i,variety) = nanmean(sig(g(:,4)==i & off));
        sig(g(:,4)==i & off) = sig(g(:,4)==i & off)./expected;
    end
    post = g(:,end)>2;
    anovabar(sig(post)-1,g(post,end)-2)
    q = [nanmean(sig(g(:,end)==3)) nanmean(sig(g(:,end)==4))];
    p = ranksum(sig(g(:,end)==3),sig(g(:,end)==4));
    text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    title(['p_d=' num2str(p)]);
    p = [signrank(sig(g(:,end)==3)-1,0,'tail','right') signrank(sig(g(:,end)==4)-1,0,'tail','right')];
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['current,p=' num2str(p(1))],['prev,p=' num2str(p(2))]},'box','off');
    set(gca,'yticklabel',get(gca,'ytick')+1);
    set(gca,'yticklabel',get(gca,'ytick')+1);
    ylabel('proportion replay (post/baseline)');
    ylim([0.9 1.2]-1);
    drawnow


    subplot(2,6,6);
    ns = Accumulate(g0(:,end),g0(:,2)<0.05 & g0(:,3)~=0);
    q = Accumulate(g0(:,end),g0(:,2)<0.05 & g0(:,3)>0);
    z = zBinomialComparison(q(3),ns(3),0.5);
    z(2) = zBinomialComparison(q(4),ns(4),0.5);
    %     prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
    p = z2p(zBinomialComparison(q(4),ns(4),q(3),ns(3)));
    q = q./ns; q = q(3:4);%./q(1:2);
    bar(q);
    text(1,(q(1))/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    text(2,(q(2))/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    title(['p_d=' num2str(p)]);
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['current,p=' num2str(z2p(z(1)))],['prev,p=' num2str(z2p(z(2)))]},'box','off');
    ylabel('proportion forward replay (in post)');
    drawnow
    ylim([0 1]);

    % normalize current vs previous (rather than with preplay)

    subplot(2,6,12);
    sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan;
    on = ismember(g(:,end),[1 3]);   off = ismember(g(:,end),[2 4]);
    pre = ismember(g(:,end),[1 2]);   post = ismember(g(:,end),[3 4]);
    for i=1:max(seshID),
        ok = g(:,end)==2 & g(:,4)==i; expected = nanmean(sig(ok)); sig(g(:,4)==i & pre) = sig(g(:,4)==i & pre)./expected;
        ok = g(:,end)==4 & g(:,4)==i; expected = nanmean(sig(ok)); sig(g(:,4)==i & post) = sig(g(:,4)==i & post)./expected;
    end

    anovabar(sig(on)-1,g(on,end)-post(on));
    ylim([-1 1]);
    q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==3))];
    p = ranksum(sig(g(:,end)==1),sig(g(:,end)==3));
    text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    title(['p_d=' num2str(p)]);
    p = z2p(z); p(z<0) = 1;
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['presleep,p=' num2str(p(1))],['postsleep,p=' num2str(p(2))]},'box','off');
    set(gca,'yticklabel',get(gca,'ytick')+1);
    set(gca,'yticklabel',get(gca,'ytick')+1);
    ylabel('proportion replay (current/previous)');
    drawnow
end

%%
clear maps;
xlimit = 10;
for i=1:max(seshID)
    figure(i);
    clf
    o = abs(g0(:,3))<30 & g(:,4)==i; o1 = o & g(:,end)==1; o2 = o & g(:,end)==3;
    subplot(2,3,1);
    [h1,x,y] = DensityMap(g0(o1,3),g0(o1,1),'nBins',[600 50],'smooth',1,'show','off');hh1 = h1./sum(h1(:));
    PlotColorMap(hh1,'x',x,'y',y); clim([0 1]*0.002);
    PlotHVLines(0,'v','k'); xlim([-1 1]*xlimit);
    subplot(2,3,2);
    h2 = DensityMap(g0(o2,3),g0(o2,1),'nBins',[600 50],'smooth',1,'show','off');hh2 = h2./sum(h2(:));
    PlotColorMap(hh2,'x',x,'y',y); clim([0 1]*0.002);
    PlotHVLines(0,'v','k'); xlim([-1 1]*xlimit);

    title(X{i,end});

    subplot(2,3,3);
    PlotColorMap(hh2-hh1,'x',x,'y',y);
    PlotHVLines(0,'v','k'); xlim([-1 1]*xlimit);
    clim([-1 1]*0.001);
    maps{i,1} = hh1; maps{i,2} = hh2; maps{i,3} = hh2-hh1;

    o = abs(g0(:,3))<30 & g(:,4)==i; o1 = o & g(:,end)==2; o2 = o & g(:,end)==4;
    subplot(2,3,4);
    [h1,x,y] = DensityMap(g0(o1,3),g0(o1,1),'nBins',[600 50],'smooth',1,'show','off');hh1 = h1./sum(h1(:));
    PlotColorMap(hh1,'x',x,'y',y); clim([0 1]*0.002);
    PlotHVLines(0,'v','k'); xlim([-1 1]*xlimit);
    subplot(2,3,5);
    h2 = DensityMap(g0(o2,3),g0(o2,1),'nBins',[600 50],'smooth',1,'show','off');hh2 = h2./sum(h2(:));
    PlotColorMap(hh2,'x',x,'y',y); clim([0 1]*0.002);
    PlotHVLines(0,'v','k'); xlim([-1 1]*10);
    subplot(2,3,6);
    PlotColorMap(hh2-hh1,'x',x,'y',y);
    PlotHVLines(0,'v','k'); xlim([-1 1]*xlimit);
    clim([-1 1]*0.001);

    maps{i,4} = hh1; maps{i,5} = hh2; maps{i,6} = hh2-hh1;
end
if conditionX==1, maps1 = maps; else maps0 = maps; end

close all
figure(20);
for i=1:6,   subplot(2,3,i); PlotColorMap(nanmean(cat(3,maps{:,i}),3),'x',x,'y',y); PlotHVLines(0,'v','k'); xlim([-1 1]*xlimit);
    if rem(i,3)==0, clim([-1 1]*0.001); else clim([0 1]*0.002); end
end

% SaveFig(fullfile('M:\home\raly\results\OML\Cheese-',[sessionID '-summary-results-excluding-678']))

%%
smooth = 0;
for conditionX = 0:1
    X = X0(cell2mat(X0(:,9))==conditionX,:);
    clear zTraj zPool
    clear savedMaps savedMapsSh savedMapsZ
    for ses=1:size(X,1)
        for q=1:2
            pre = X{ses,5}; post = X{ses,6};
            outputs = X{ses,q+2}; % outputsCurrent;
            % outputs = X{i,4}; % outputsPrev;
            % outputs = outputsPrev;
            for i=1:size(outputs,1),for j=1:size(outputs,2), if isempty(outputs{i,j}), outputs{i,j} = nan(size(outputs{1,j})); end; end; end
            outputs1 = outputs(pre,:);

            % figure
            score = cell2mat(outputs1(:,1));
            scoreSh = cell2mat(outputs1(:,5)')';
            slope = cell2mat(outputs1(:,4))-cell2mat(outputs1(:,3));
            slopeSh = cell2mat(outputs1(:,7)')'-cell2mat(outputs1(:,6)')';
            nShuffles = size(scoreSh,2);
            jump = cell2mat(outputs1(:,10));
            jumpSh = cell2mat(outputs1(:,11)')';
            c = cell2mat(outputs1(:,8));
            cSh = cell2mat(outputs1(:,9)')';
            silva = c>0.4 & jump<0.4*100;
            silvaSh = cSh>0.4 & jumpSh<0.4*100;
            zTraj(ses,1+(q-1)*2) = (sum(silva)-mean(sum(silvaSh)))./std(sum(silvaSh));
            %     [sum(silva) sum(~isnan(c))]

            varX = jump; varXsh = jumpSh; binsX = Bins(0,100,5);
            varY = c; varYsh = cSh; binsY = Bins(0,1,0.0333);
            varY = score; varYsh = scoreSh; binsY = Bins(0,1,1/100);
            varX = slope; varXsh = slopeSh; binsX = Bins(-200,200,400/100);

            x = mean(binsX,2)/100; y = mean(binsY,2);
            [~,wy] = InIntervals(varY,binsY); [~,wySh] = InIntervals(varYsh(:),binsY);
            [~,wx] = InIntervals(varX,binsX); [~,wxSh] = InIntervals(varXsh(:),binsX);
            wySh = reshape(wySh,size(varYsh)); wxSh = reshape(wxSh,size(varXsh));
            ok = ~isnan(varX) & ~isnan(varY) & wy>0 & wx>0;
            h = Accumulate([wy(ok) wx(ok)],1,'size',[size(binsY,1) size(binsX,1)]);
            h = Smooth(h,smooth);
            hSh = nan([size(h) nShuffles]);
            for i=1:nShuffles,
                hSh(:,:,i) = Smooth(Accumulate([wySh(ok,i) wxSh(ok,i)],1,'size',[size(binsY,1) size(binsX,1)]),smooth);
            end
            z = (h-nanmean(hSh,3))./nanstd(hSh,[],3);

            savedMaps{ses,q+0} = h;
            savedMapsSh{ses,q+0} = hSh;
            savedMapsZ{ses,q+0} = z;

            %
            %             subplot(2,4,1);
            %             PlotColorMap(h,'x',x,'y',y);
            %             subplot(2,4,2);
            %             PlotColorMap(nanmean(hSh,3),'x',x,'y',y);
            %             subplot(2,4,3);
            %             PlotColorMap(z,'x',x,'y',y);
            %             % set(gca,'ydir','reverse');
            % %             clim([0 z2p(0.01)]);
            %             subplot(2,4,4);
            %             PlotColorMap(h-nanmean(hSh,3),'x',x,'y',y);
            %             hPre = h/nanmean(h(:)); dPre = (h-nanmean(hSh,3))/nanmean(h(:));

            % post
            outputs1 = outputs(post,:);

            % figure
            score = cell2mat(outputs1(:,1));
            scoreSh = cell2mat(outputs1(:,5)')';
            slope = cell2mat(outputs1(:,4))-cell2mat(outputs1(:,3));
            slopeSh = cell2mat(outputs1(:,7)')'-cell2mat(outputs1(:,6)')';
            nShuffles = size(scoreSh,2);
            jump = cell2mat(outputs1(:,10));
            jumpSh = cell2mat(outputs1(:,11)')';
            c = cell2mat(outputs1(:,8));
            cSh = cell2mat(outputs1(:,9)')';
            silva = c>0.4 & jump<0.4*100;
            silvaSh = cSh>0.4 & jumpSh<0.4*100;
            zTraj(ses,2+(q-1)*2) = (sum(silva)-mean(sum(silvaSh)))./std(sum(silvaSh));
            %     [sum(silva) sum(~isnan(c))]
            
            varX = jump; varXsh = jumpSh;
            varY = c; varYsh = cSh; 
            varY = score; varYsh = scoreSh; 
            varX = slope; varXsh = slopeSh; 

            [~,wy] = InIntervals(varY,binsY); [~,wySh] = InIntervals(varYsh(:),binsY);
            [~,wx] = InIntervals(varX,binsX); [~,wxSh] = InIntervals(varXsh(:),binsX);
            wySh = reshape(wySh,size(varYsh)); wxSh = reshape(wxSh,size(varXsh));
            ok = ~isnan(varX) & ~isnan(varY);
            h = Accumulate([wy(ok) wx(ok)],1,'size',[size(binsY,1) size(binsX,1)]);
            h = Smooth(h,smooth);
            hSh = nan([size(h) nShuffles]);
            for i=1:nShuffles,
                hSh(:,:,i) = Smooth(Accumulate([wySh(ok,i) wxSh(ok,i)],1,'size',[size(binsY,1) size(binsX,1)]),smooth);
            end
            z = (h-nanmean(hSh,3))./nanstd(hSh,[],3);


            savedMaps{ses,q+2} = h;
            savedMapsSh{ses,q+2} = nanmean(hSh,3);
            savedMapsZ{ses,q+2} = z;

        end
        ses
    end
    if conditionX==0
        savedOFF = {savedMaps,savedMapsSh,savedMapsZ};
    else
        savedON = {savedMaps,savedMapsSh,savedMapsZ};
    end
end

%%

for conditionX = 0:1
    if conditionX==0,
        [savedMaps,savedMapsSh,savedMapsZ] = savedOFF{:};
    else
        [savedMaps,savedMapsSh,savedMapsZ] = savedON{:};
    end
    for j=1:2,
        for q=1:2,
            if conditionX==0, name = 'off'; else name = 'stim'; end
            if j==1, name = [name ', pre']; else name = [name ', post']; end
            if q==1, name = [name ', current']; else name = [name ', previous']; end
            subplot(2,4,j+(q-1)*4 + 2*conditionX);
            this = nanmean(cat(3,savedMapsZ{:,q+(j-1)*2}),3); this(isnan(this)) = 0; this(this==Inf) = 100;
            PlotColorMap(Smooth(this,2),'x',x,'y',y);grid on; clims([0 10]);
            xlim([-1 1]);
            title(name)
        end
    end
end


%%
% pooled

for q=1:2
    pre = cat(1,X{:,5}); post = cat(1,X{:,6});
    outputs = cat(1,X{:,q+2}); % outputsCurrent;
    % outputs = X{i,4}; % outputsPrev;
    % outputs = outputsPrev;
    for i=1:size(outputs,1),for j=1:size(outputs,2), if isempty(outputs{i,j}), outputs{i,j} = nan(size(outputs{1,j})); end; end; end
    outputs1 = outputs(~post,:);
    jump = cell2mat(outputs1(:,10));
    jumpSh = cell2mat(outputs1(:,11)')';
    c = cell2mat(outputs1(:,8));
    cSh = cell2mat(outputs1(:,9)')';
    silva = c>0.4 & jump<0.4*100;
    silvaSh = cSh>0.4 & jumpSh<0.4*100;
    zPool(1,1+(q-1)*2) = (sum(silva)-mean(sum(silvaSh)))./std(sum(silvaSh));
    Portion(silva)
    % post
    outputs1 = outputs(post,:);

    % figure
    jump = cell2mat(outputs1(:,10));
    jumpSh = cell2mat(outputs1(:,11)')';
    c = cell2mat(outputs1(:,8));
    cSh = cell2mat(outputs1(:,9)')';
    silva = c>0.4 & jump<0.4*100;
    silvaSh = cSh>0.4 & jumpSh<0.4*100;
    zPool(1,2+(q-1)*2) = (sum(silva)-mean(sum(silvaSh)))./std(sum(silvaSh));
    Portion(silva)
end

% z1 = [zTraj; zPool];


%%





