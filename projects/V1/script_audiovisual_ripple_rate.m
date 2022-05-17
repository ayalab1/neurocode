% script_audiovisual_ripple_rate

X = BatchReturn('AudioVisual.batch');

for i=1:size(X,1)
    basepath = X{i,1};

    ripples = getStruct(basepath,'ripples');
    MergePoints = getStruct(basepath,'MergePoints');
    EMG = getStruct(basepath,'EMG');
    immobility = EMG.timestamps(FindInterval(EMG.data<0.6)); immobility(diff(immobility,[],2)<1,:) = [];
    postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
    sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));

    r = ripples.timestamps;
    r = Restrict(r,immobility);
    pre = InIntervals(r(:,1),sleep(1,:));
    post = InIntervals(r(:,1),sleep(2,:));
    rippleRate(i,:) = [sum(pre)./dur(Restrict(immobility,sleep(1,:))) sum(post)./dur(Restrict(immobility,sleep(2,:)))];
end

rippleRate(any(isnan(rippleRate),2),:) = [];
clf
anovaplot(rippleRate,[],'parametric',false,'alpha',[0 0.05]);
hold all
plot(meshgrid(1:size(rippleRate,2),1:size(rippleRate,1))'+randn(size(rippleRate'))/20,rippleRate','k.-','markersize',20,'linewidth',1.5)
xlim([-1 3]);
set(gca,'xtick',[1 2],'xticklabel',{'pre-task sleep','post-task sleep'});
ylabel('Ripple rate (Hz)');
set(gca,'fontsize',12,'tickdir','out','box','off');


a = Insets('left',[0.2 0.4],'margin',0.1);
q = bootstrp(100000,@(x) mean(diff(x,[],2)),rippleRate);
[h,ht] = hist(q,linspace(-0.2,0.5,1000)); 
ci = quantile(q,[0.025 0.975]);
colors = get(gca,'ColorOrder');
% plot(Smooth(h,10),ht,'color',colors(1,:)); hold on
handle = fill(Smooth(h,10),ht,colors(1,:));
set(handle,'facealpha',0.5,'EdgeColor','none');
hold on
% plot(-Smooth(h,10),ht,'color',colors(1,:));
plot([0 0],ci,'k','linewidth',3);
plot(0,median(q),'k.','markersize',30);
PlotHVLines(0,'h','k--');
set(gca,'box','off','xtick',[]);
xlim([-500 1000]);
ylabel('post-pre difference (Hz)');


SaveFig('M:\home\raly\results\V1\Pooled\Audiovisual_training_rippleRate');