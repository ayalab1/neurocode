batch = StartBatch(@BatchCanReplay,'OMLproject.batch');

X = get(batch,'UserData');
%%


eon = cat(2,X{:,1});
eoff = cat(2,X{:,2});

sessionID = 'pooled';
on = squeeze(nanmean(reshape(eon,200,6,[]),3));
off = squeeze(nanmean(reshape(eoff,200,6,[]),3));

%%

clf
subplot(1,4,1);
average = on;
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
PlotColorMap(repmat(average,1,2));
hold on; plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep([' stim ON'],'_','-'));

subplot(1,4,2);
average = off;
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
PlotColorMap(repmat(average,1,2));
hold on;
plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep([' stim OFF'],'_','-'));

subplot(1,2,2);
Qok = false(200,6); Qok(75:99,2:3) = true; Qok(101:125,4:5) = true;
Qcontrol = false(200,6); Qcontrol(75:99,4:5) = true; Qcontrol(101:125,2:3) = true;
this = reshape(eon,200,6,[]); this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
this = reshape(eon,200,6,[]); this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
onScore = score;
this = reshape(eoff,200,6,[]); this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
this = reshape(eoff,200,6,[]); this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
offScore = score;
g = Group(onScore,offScore); g(isnan(g(:,2)),:) = [];
gg =[-diff(g(:,1:2)./sum(g(:,1:2),2),[],2),g(:,3)];
anovaplot(-diff(g(:,1:2)./sum(g(:,1:2),2),[],2),g(:,3),'parametric','off');
p = [signrank(gg(gg(:,2)==1)) signrank(gg(gg(:,2)==2))];
set(gca,'xtick',1:2,'xticklabel',{['ON, p=' num2str(p(1))],['OFF, p=' num2str(p(2))]});
p = ranksum(gg(gg(:,2)==1),gg(gg(:,2)==2));
title(['Quadrant scores: ranksum p=' num2str(p)]);

%%

outputsOn = cat(1,X{:,3}); outputsOff = cat(1,X{:,4}); outputsBoth = cat(1,X{:,5});
pre = cat(1,X{:,6}); post = cat(1,X{:,7}); durations = cat(1,X{:,8});
% Need to normalize relative to pre-sleep separately for each session!!!

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
        these{condition,k} = [scores pValues slopes];
    end
end

g = Group(these{:});
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));

subplot(3,2,1);
ok = g(:,end)>2;
anovabar(gg(ok,1),g(ok,end)-2,'parametric',false)
title([strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');
% subplot(2,2,2);
% q = Accumulate(g(:,end),g(:,2)<0.05)./Accumulate(g(:,end),1); 
% q = Accumulate(g(:,end),g(:,3)>0)./Accumulate(g(:,end),1); 
% q = q(3:4)./q(1:2); bar(q-1);
% title(['all bursts: ' num2str(round(q(1)*1000)/1000) ', ' num2str(round(1000*q(2))/1000)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
% set(gca,'yticklabel',get(gca,'ytick')+1);
% ylabel('proportion replay (post/baseline)');
% 

subplot(3,4,3);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');

subplot(3,4,4);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');

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
        these{condition,k} = [scores pValues slopes];
    end
end

g = Group(these{:});
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));

subplot(3,2,3);
ok = g(:,end)>2;
anovabar(gg(ok,1),g(ok,end)-2,'parametric',false)
title(['better fit (on/off): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');

subplot(3,4,7);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');


subplot(3,4,8);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');

% cla
% anovabar(g(:,3)>0,g(:,end));


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
        these{condition,k} = [scores pValues slopes];
    end
end

g = Group(these{:});
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));

subplot(3,2,5);
ok = g(:,end)>2;
anovabar(gg(ok,1),g(ok,end)-2,'parametric',false)
title(['decoded together (sorted posthoc): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');

subplot(3,4,11);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');


subplot(3,4,12);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');


subplot(3,2,1); y = ylim; subplot(3,2,3); y = [min([min(ylim) y(1)]) max([max(ylim) y(2)])]; ylim(y); subplot(3,2,1); ylim(y);
subplot(3,4,3); y = ylim; 
for k=[4 7 8 11 12], subplot(3,4,k); y = [min([min(ylim) y(1)]) max([max(ylim) y(2)])]; end
for k=[3 4 7 8 11 12], subplot(3,4,k); ylim(y); set(gca,'yticklabel',get(gca,'ytick')+1); end


