batch = StartBatch(@BatchV1ProbeTrials,'AudioVisual.batch');
X = get(batch,'UserData');

%%
clf
g = cell2mat(X(:,1));
hs = cell2mat(X(:,2));
ht = linspace(-1,1,201)*0.5;
handle0 = semplot(ht,hs(hs(:,end)==3,1:end-1),'k');
handle1 = semplot(ht,hs(hs(:,end)==1,1:end-1),'r');
handlep = semplot(ht,hs(hs(:,end)==2,1:end-1),'b');
PlotHVLines([0.09 0.14],'v','k--','linewidth',2);
legend([handle0;handle1;handlep],'baseline','first 10 responses (non-probe)','worst putative probe trials');
legend('box','off');
set(gca,'fontsize',12,'tickdir','out','box','off');
xlabel('time from visual cue presentation (s)');
ylabel('V1 multiunit response (Hz)');

a = Insets('top',[0.2 0.4],'margin',0.1);
anovaplot(g(:,1),g(:,2),'alpha',[0 0.05]);
ylim([-20 300]);
set(a,'xtick',1:3,'xticklabel',{'first10','worst','ctrl'});
set(a,'fontsize',12,'tickdir','out','box','off');

SaveFig('M:\home\raly\results\V1\Pooled\PutativeProbeTrialResponse');