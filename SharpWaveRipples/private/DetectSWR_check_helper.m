% DetectSWR_check_helper - script that can be manually called from DetectSWR

% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
    [clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[8 5],'manual',true);
    bad = bad | InIntervals(t,badIntervals);
end
if false % optional steps (not recommended for the general case)
    try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
        % if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
        % reasonable to you!
        smoothed = Smooth(abs(swDiff),1250*10);
        noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
        bad = bad | InIntervals(t,noisyPeriods);
    end
    try % optionally, remove ripples in which the sharp wave channel was positive
        % Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
        sw = interp1(tl,lfpLow(:,2),t);
        %                 bad = bad | sw>0;
    end
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"

figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));


scores = nan(size(ripPowerAll,1),1);
done = false;
while ~done
    figure(1);
    colors = Bright(1000);
    [x,y,button] = ginput(1);
    [~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));

    t = featureTs/1250;
    tl = (1:length(lfp))'/1250;
    figure(2);
    clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
    subplot(3,1,1);
    plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
    hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
    legend('ripple channel','SW channel','noise channel');
    subplot(3,1,2);
    plot(tl(in) - t(i),double(swDiff(in)));
    hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
    subplot(3,1,3);
    plot(tl(in) - t(i),double(ripPower0(in)));
    hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
    drawnow

    %     x = input( prompt )
    commandwindow % select the command window
    answer = input('Score as good(1)/bad(0) or skip (hit Enter), exit (type ''done''), or enter debug mode (type ''keyboard''): ','s');
    if strcmp(lower(answer),'keyboard')
        keyboard
    end
    if strcmp(lower(answer),'done')
        done = true;
    end
    score = str2double(answer);
    scores(i) = score; % save scores to optionally save

    figure(1);
    xlims = xlim; ylims = ylim; clf
    % extrapolate the score of each ripple based on the score of the nearest scored ripple
    scored = find(~isnan(scores));
    distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
    [~,nearestNeighbour] = min(distances,[],2);
    estimated = scores(scored(nearestNeighbour));
    selected = estimated>0.1;
    idx1 = selected & ~bad; % final ripples
    idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)

    scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
    xlim(xlims); ylim(ylims);
    hold on;
    scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
    set(get(colorbar,'YLabel'),'String','Estimated ripple score');
    xlabel('Sharp wave depth'); ylabel('Ripple power');

    title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});

    if rem(i,10)==0
        display('saving...');
        save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
    end
end

if false % optional steps (not recommended for the general case)
    % Otherwise, you can just draw a polyglon:
    selected = UIInPolygon(swDiffAll,ripPowerAll); % draw polygon to encircle the points you believe are ripples
    idx1 = selected & ~bad; % final ripples
    idx2 = ~selected & ~bad;
end
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
