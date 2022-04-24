batch = StartBatch(@BatchCanReplay,'OMLproject.batch');

X = get(batch,'UserData');
%%


% eon = cat(2,X{:,1});
% eoff = cat(2,X{:,2});

sessionID = 'pooled';
on = squeeze(nanmean(reshape(eon,200,6,[]),3));
off = squeeze(nanmean(reshape(eoff,200,6,[]),3));

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

% %%
% clf
% Qok = false(200,6); Qok(75:99,2:3) = true; Qok(101:125,4:5) = true;
% Qcontrol = false(200,6); Qcontrol(75:99,4:5) = true; Qcontrol(101:125,2:3) = true;
% this = reshape(eon,200,6,[]); this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
% this = reshape(eon,200,6,[]); this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
% onScore = score;
% this = reshape(eoff,200,6,[]); this(~repmat(Qok,1,1,size(this,3))) = nan; score = reshape(nanmean(nanmean(this,1),2),[],1);
% this = reshape(eoff,200,6,[]); this(~repmat(Qcontrol,1,1,size(this,3))) = nan; score(:,2) = reshape(nanmean(nanmean(this,1),2),[],1);
% offScore = score;
% g = Group(onScore,offScore); g(isnan(g(:,2)),:) = [];
% gg =[-diff(g(:,1:2)./sum(g(:,1:2),2),[],2),g(:,3)];
% anovaplot(-diff(g(:,1:2)./sum(g(:,1:2),2),[],2),g(:,3),'parametric','off');
% set(gca,'xtick',1:2,'xticklabel',{'ON','OFF'});

%%

outputsOn,outputsOff,



