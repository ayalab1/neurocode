batch = StartBatch(@BatchCanReplay,'OMLproject.batch');

X0 = get(batch,'UserData');
%%
figure(1); 
sessionID = 'O';
X = X0; 
ok = cellfun(@(x) ~isempty(strfind(x,sessionID)),X0(:,end));

X(~ok,:) = [];  %X(14:16,:) = [];
eon = cat(2,X{:,1});
eoff = cat(2,X{:,2});

mon = squeeze(nanmean(reshape(eon,200,6,[]),3));
moff = squeeze(nanmean(reshape(eoff,200,6,[]),3));

% figure('name','18');
clf
subplot(2,6,1);
average = mon;
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
PlotColorMap(repmat(average,1,2));
hold on; plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep(['Theta: stim ON'],'_','-'));

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
title(strrep(['Theta: stim OFF'],'_','-'));
clims
drawnow

subplot(2,6,7); cla
sessionNames = cellfun(@(x) x(9:end),X(:,end),'UniformOutput',0);
text(0,0.5,sessionNames);
set(gca,'xtick',[],'ytick',[])
title('sessions included');

subplot(2,6,8);
Qok = false(200,6); Qok(75:99,2:3) = true; Qok(101:125,4:5) = true;
Qcontrol = false(200,6); Qcontrol(75:99,4:5) = true; Qcontrol(101:125,2:3) = true;
this = reshape(eon,200,6,[]); this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
this = reshape(eon,200,6,[]); this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
onScore = score;
this = reshape(eoff,200,6,[]); this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
this = reshape(eoff,200,6,[]); this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
offScore = score;
g0 = Group(onScore,offScore); g0(isnan(g0(:,2)),:) = [];
g =[-diff(g0(:,1:2)./sum(g0(:,1:2),2),[],2),g0(:,3)];
anovaplot(-diff(g0(:,1:2)./sum(g0(:,1:2),2),[],2),g0(:,3),'parametric','off');
p = [signrank(g(g(:,2)==1)) signrank(g(g(:,2)==2))];
set(gca,'xtick',1:2,'xticklabel',{['ON, p=' num2str(p(1))],['OFF, p=' num2str(p(2))]});
p = ranksum(g(g(:,2)==1),g(g(:,2)==2));
title(['Theta: quadrant scores: ranksum p=' num2str(p)]);

% Replay
outputsOn = cat(1,X{:,3}); outputsOff = cat(1,X{:,4}); outputsBoth = cat(1,X{:,5});
pre = cat(1,X{:,6}); post = cat(1,X{:,7}); durations = cat(1,X{:,8});
seshID = repelem((1:size(X,1))',cellfun(@length,X(:,6)),1);

for variety = 1:3,
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
        case 3 % decode together:
            outputs = outputsBoth;
            reOn = false(size(outputs,1),1);
            reOff = false(size(outputs,1),1);
            for i=1:size(outputsOn,1), try reOn(i,1) = outputs{i,3}<=size(average,1)/2; end; end
            for i=1:size(outputsOn,1), try reOff(i,1) = outputs{i,3}>size(average,1)/2; end; end
            for condition = 1:2
                if condition==1
                    bad = cellfun(@isempty,outputs(:,1));
                    bad = bad | ~(reOn);
                else
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

    subplot(3,3,(variety-1)*3+2);
    ok = g0(:,end)>2;
    anovabar(g(ok,1),g0(ok,end)-2,'parametric',false)
    if variety==1
        title(['Replay: ' strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
    elseif variety==2
        title(['Replay: better fit (on/off): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
    else
        title(['Replay: decoded together (sorted posthoc): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
    end
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
    ylabel('score relative to baseline');

    drawnow

    subplot(3,6,(variety-1)*6+5);
    % normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
    sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; on = ismember(g(:,end),[1 3]);   off = ismember(g(:,end),[2 4]);
    for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); sig(g(:,4)==i & on) = sig(g(:,4)==i & on)./expected; 
        preOn(i,variety) = expected; postOn(i,variety) = nanmean(sig(g(:,end)==3 & g(:,4)==i));
        ok = g(:,end)==2 & g(:,4)==i; expected = nanmean(sig(ok)); sig(g(:,4)==i & off) = sig(g(:,4)==i & off)./expected; 
        
    end
    post = g(:,end)>2; 
    anovabar(sig(post)-1,g(post,end)-2)
    q = [nanmean(sig(g(:,end)==3)) nanmean(sig(g(:,end)==4))];
    p = ranksum(sig(g(:,end)==3),sig(g(:,end)==4));
    text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
    title(['p_d=' num2str(p)]);
    p = [signrank(sig(g(:,end)==3)-1) signrank(sig(g(:,end)==4)-1)];
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(p(1))],['off,p=' num2str(p(2))]},'box','off');
    set(gca,'yticklabel',get(gca,'ytick')+1);
    set(gca,'yticklabel',get(gca,'ytick')+1);
    ylabel('proportion replay (post/baseline)');
    drawnow

    subplot(3,6,(variety-1)*6+6);
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
    set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
    set(gca,'yticklabel',get(gca,'ytick')+1);
    ylabel('proportion forward replay (in post)');
    drawnow
end

% SaveFig(fullfile('M:\home\raly\results\OML\',[sessionID '-summary-results-excluding-678']))

%%



% 
% g0 = Group(these{:});
% g = g0; g(ismember(g(:,end),[1 3]),1) = (g(ismember(g(:,end),[1 3]),1) - nanmedian(g(g(:,end)==1,1)))./diff(quantile(g(g(:,end)==1,1),[0.25 0.75]));
% g(ismember(g(:,end),[2 4]),1) = (g(ismember(g(:,end),[2 4]),1) - nanmedian(g(g(:,end)==2,1)))./diff(quantile(g(g(:,end)==2,1),[0.25 0.75]));
% 
% subplot(3,2,3);
% ok = g0(:,end)>2;
% anovabar(g(ok,1),g0(ok,end)-2,'parametric',false)
% title(['better fit (on/off): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
% ylabel('score relative to baseline');
% 
% subplot(3,4,7);
% ns = Accumulate(g0(:,end),1);
% q = Accumulate(g0(:,end),g0(:,2)<0.05);
% z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
% z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
% prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
% p = z2p(zBinomialComparison(q(4),ns(4),prediction));
% q = q./ns; q = q(3:4)./q(1:2); 
% bar(q-1);
% text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% title(['p_d=' num2str(p)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
% set(gca,'yticklabel',get(gca,'ytick')+1);
% ylabel('proportion replay (post/baseline)');
% 
% 
% subplot(3,4,8);
% ns = Accumulate(g0(:,end),1);
% q = Accumulate(g0(:,end),g0(:,2)<0.05 & g0(:,3)>0);
% z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
% z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
% prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
% p = z2p(zBinomialComparison(q(4),ns(4),prediction));
% q = q./ns; q = q(3:4)./q(1:2); 
% bar(q-1);
% text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% title(['p_d=' num2str(p)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
% set(gca,'yticklabel',get(gca,'ytick')+1);
% ylabel('proportion forward replay (post/baseline)');
% 
% % cla
% anovabar(g(:,3)>0,g(:,end));


% BOTH
% g0 = Group(these{:});
% g = g0; g(ismember(g(:,end),[1 3]),1) = (g(ismember(g(:,end),[1 3]),1) - nanmedian(g(g(:,end)==1,1)))./diff(quantile(g(g(:,end)==1,1),[0.25 0.75]));
% g(ismember(g(:,end),[2 4]),1) = (g(ismember(g(:,end),[2 4]),1) - nanmedian(g(g(:,end)==2,1)))./diff(quantile(g(g(:,end)==2,1),[0.25 0.75]));
% 
% subplot(3,2,5);
% ok = g0(:,end)>2;
% anovabar(g(ok,1),g0(ok,end)-2,'parametric',false)
% title(['decoded together (sorted posthoc): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
% ylabel('score relative to baseline');
% 
% subplot(3,4,11);
% ns = Accumulate(g0(:,end),1);
% q = Accumulate(g0(:,end),g0(:,2)<0.05);
% z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
% z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
% prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
% p = z2p(zBinomialComparison(q(4),ns(4),prediction));
% q = q./ns; q = q(3:4)./q(1:2); 
% bar(q-1);
% text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% title(['p_d=' num2str(p)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
% set(gca,'yticklabel',get(gca,'ytick')+1);
% ylabel('proportion replay (post/baseline)');
% 
% 
% subplot(3,4,12);
% ns = Accumulate(g0(:,end),1);
% q = Accumulate(g0(:,end),g0(:,2)<0.05 & g0(:,3)>0);
% z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
% z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
% prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
% p = z2p(zBinomialComparison(q(4),ns(4),prediction));
% q = q./ns; q = q(3:4)./q(1:2); 
% bar(q-1);
% text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
% title(['p_d=' num2str(p)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
% set(gca,'yticklabel',get(gca,'ytick')+1);
% ylabel('proportion forward replay (post/baseline)');
% 
% 
% subplot(3,2,1); y = ylim; subplot(3,2,3); y = [min([min(ylim) y(1)]) max([max(ylim) y(2)])]; ylim(y); subplot(3,2,1); ylim(y);
% subplot(3,4,3); y = ylim; 
% for k=[4 7 8 11 12], subplot(3,4,k); y = [min([min(ylim) y(1)]) max([max(ylim) y(2)])]; end
% for k=[3 4 7 8 11 12], subplot(3,4,k); ylim(y); set(gca,'yticklabel',get(gca,'ytick')+1); end


