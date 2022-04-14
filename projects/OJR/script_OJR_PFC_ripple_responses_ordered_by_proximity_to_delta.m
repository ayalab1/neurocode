batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
Xs = get(batchS,'UserData'); % short
Xl = get(batchL,'UserData'); % long
Xo = get(batchO,'UserData'); % opto

Xs(cellfun(@isempty,Xs(:,1)),:) = []; Xl(cellfun(@isempty,Xl(:,1)),:) = []; Xo(cellfun(@isempty,Xo(:,1)),:) = [];
XX = {Xs,Xl,Xo};
write = false;

%%
% figure;
clf
for condition=1:3
    X = XX{condition};
    for i=1:size(X,1)

        [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:}; 
        ripples = Restrict(ripples,sws);
        if isempty(pfc), continue; end
        try
            upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(~InIntervals(diff(upstate,[],2),[1 4]),:) = [];
        catch
            continue;
        end
        rt = RelativeTime(ripples(:,1),upstate);

        [h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
        [hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
 
            subplot(3,5,(condition-1)*5 + i)
        
        rt = ripples(:,1) - deltas(FindClosest(deltas(:,2),ripples(:,1)),2);
        ok = InIntervals(rt,[-1 1]*0.5);
        [~,m] = max(Shrink(sortby(hd(ok,:),rt(ok)),floor(sum(ok)/100),1),[],2);
        PlotColorMap(Smooth(Shrink(sortby(h(ok,:),rt(ok)),floor(sum(ok)/100),1),2),'x',ht);
        hold all
        plot(ht(m),1:max(ylim),'w','linewidth',2);
        PlotHVLines(0,'v','k--','linewidth',2);
        title(basepath);
        drawnow
    end
end

if write
    SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_proximity_to_delta']));
end

% [PETH1,t1] = PETH(ripples(:,2),deltas(:,2),'durations',[-1 1]*0.5,'nBins',201); PETH1 = PETH1/mode(diff(t1)); % in Hz
% [PETH2,t2] = PETH(spikes(:,1),deltas(:,2),'durations',[-1 1]*0.5,'nBins',201); PETH2 = PETH2/mode(diff(t2)); % in Hz
% [joint, expected, difference] = JointPETH(PETH1,PETH2,2);
% % figure;
% subplot(2,2,1); PlotColorMap(joint,'x',t2,'y',t1);
% ylabel('ripple rate'); xlabel('spike rate');
% PlotHVLines(0,'b','k--','linewidth',2);
% subplot(2,2,2); PlotColorMap(expected,'x',t2,'y',t1);
% subplot(2,2,3); PlotColorMap(difference,'x',t2,'y',t1);
% PlotHVLines(0,'b','k--','linewidth',2);

%%

% figure;
clf
nBins = 100; smooth = 1;
for condition=1:3
    X = XX{condition};
    for i=1:size(X,1)

        [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
        ripples = Restrict(ripples,sws);
        if isempty(pfc) || isempty(deltas) || isempty(ripples), continue; end
        try
            upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(~InIntervals(diff(upstate,[],2),[1 5]),:) = [];
        catch
            continue;
        end
        rng(0);
        control = sort(ripples(:,1)+1+rand(size(ripples(:,1)))*9);
        [h,ht] = PETH(pfc(:,1),control,'durations',[-1 1]*0.5,'nBins',201);
        rtC = RelativeTime(control,upstate);
        responseC = mean(h(:,InIntervals(ht,[0 0.1])),2)./0.2;

        rt = RelativeTime(ripples(:,1),upstate);
        [h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
        response = mean(h(:,InIntervals(ht,[0 0.1])),2)./0.2;
        
        subplot(3,4,(condition-1)*4 + i);
        handle = semplotBin([rt-1; rt],repmat(response,2,1),nBins,'r',smooth);
        handle0 = semplotBin([rtC-1; rtC],repmat(responseC,2,1),nBins,'k',smooth);
        set(gca,'xtick',-1:1/4:1,'xticklabel',{'delta','early upstate','mid-upstate','late upstate','delta','early upstate','mid-upstate','late upstate','delta'});
        PlotHVLines(0,'v','k--','linewidth',2);
        ylabel('multiunit PFC firing rate (Hz)');

        legend([handle handle0],{'100ms following ripples','following random control'}); legend('location','southwest','box','off');
        title(basepath);
        drawnow
    end
end
if write
    SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_relative_upstate_time']));
end
