batchL = StartBatch(@BatchCCGtransmission,'OMLlinear.batch');
batchCheese = StartBatch(@BatchCCGtransmission,'OMLcheese.batch');

%%
% subplot(2,3,1);
% ok = p(:,8)==1 & p(:,5)==0 & abs(d)<1;%; & ismember(p(:,6),[2 3]);
% g = [p(ok,1:2),p(ok,6)];% g(:,1) = zscore(g(:,1));
% anovabar(g(:,1:2),g(:,end),'parametric',false);
% ylabel('p(transmission');
% title('pooled CA3 activity -> CA1 int');
% names = {'no fields','field in stim only','field in off only','both'};
% for j=1:4, names{j} = [names{j} ',n=' num2str(sum(g(:,end)==j))]; end
% set(gca,'xtick',1:4,'xticklabel',names,'XTickLabelRotation',20);
% subplot(2,3,1+3);
% g = [abs(d(ok,1)),p(ok,6)]; %g(:,1) = zscore(g(:,1));
% anovabar(g(:,1),g(:,2),'parametric',false);
% ylabel('abs(post-pre)/(post+pre)');
% title('pooled CA3 activity -> CA1 int');
% set(gca,'xtick',1:4,'xticklabel',names,'XTickLabelRotation',20);

cheese = true;
region = 1; % 1 = CA3 -> CA1; 3 = CA3 -> CA3;
if ~cheese,
    X = get(batchL,'UserData');
else
    X = get(batchCheese,'UserData');
end

% X = X(9:end,:);
pCA3 = cell2mat(X(:,4)); pAll = cell2mat(X(:,3));
p = nani(pCA3); 
% the columns are [1: pTranmission in preSWS, 2: pTransmission in postSWS, 3: ID of source, 4: ID of target, 
% 5: is target a pyramidal cell (0 no or 1 yes), 6: is the target 0: not a pyramidal cell, 1 recorded in a stim session or 2 recorded in a control session
% 7: the region of the source activity (CA3); 8: the region of the target activity (CA1 vs CA3);
% p = nani(pAll);
q=p; 
dd = diff(q(:,1:2),[],2)./sum(q(:,1:2),2);
d = diff(q(:,1:2),[],2);
if cheese, p(:,6) = p(:,6)-1; end

% figure
clf
subplot(1,3,1);
ok = p(:,8)==region & p(:,5)>0 & abs(d)<1;% & ismember(p(:,6),[2 3]);
g = [p(ok,1:2)*100,p(ok,6)];% g(:,1) = zscore(g(:,1));
anovabar(g(:,1:2),g(:,end),'parametric',true,'alpha',[0 0.05]);
ylabel('p(transmission) (%)');
legend('pre','post sws'); legend('location','northwest','box','off','fontsize',15); 
title(['pooled CA3 activity -> CA' num2str(region) ' pyr']);
set(gca,'box','off','fontsize',15);
if ~cheese,
    names = {'no fields','field in stim only','field in off only','both'};
    for j=1:4, names{j} = [names{j} ',n=' num2str(sum(g(:,end)==j))]; end
else
    names = {'stim sessions','control sessions'};
    for j=1:2, names{j} = [names{j} ',n=' num2str(sum(g(:,end)==j))]; end
end
set(gca,'xtick',1:length(names),'xticklabel',names,'XTickLabelRotation',20);
subplot(2,3,2);
g = [(d(ok,1)),p(ok,6)]; %g(:,1) = zscore(g(:,1));
anovabar(g(:,1),g(:,2),'parametric',true,'alpha',[0 0.05]);
ylabel('(post-pre)');
title(['pooled CA3 activity -> CA' num2str(region) ' pyr']);
set(gca,'xtick',1:length(names),'xticklabel',names,'XTickLabelRotation',20);
set(gca,'box','off','fontsize',15); 

subplot(2,3,3);
g = [abs(d(ok,1)),p(ok,6)]; %g(:,1) = zscore(g(:,1));
anovabar(g(:,1),g(:,2),'parametric',false,'alpha',[0 0.05]);
ylabel('abs(post-pre)');
title(['pooled CA3 activity -> CA' num2str(region) ' pyr']);
set(gca,'xtick',1:length(names),'xticklabel',names,'XTickLabelRotation',20);
set(gca,'box','off','fontsize',15); 

subplot(2,3,5);
g = [(dd(ok,1)),p(ok,6)]; %g(:,1) = zscore(g(:,1));
anovabar(g(:,1),g(:,2),'parametric',true,'alpha',[0 0.05]);
ylabel('(post-pre)/(post+pre)');
title(['pooled CA3 activity -> CA' num2str(region) ' pyr']);
set(gca,'xtick',1:length(names),'xticklabel',names,'XTickLabelRotation',20);
set(gca,'box','off','fontsize',15); 

subplot(2,3,6);
g = [abs(dd(ok,1)),p(ok,6)]; %g(:,1) = zscore(g(:,1));
anovabar(g(:,1),g(:,2),'parametric',false,'alpha',[0 0.05]);
ylabel('abs(post-pre)/(post+pre)');
title(['pooled CA3 activity -> CA' num2str(region) ' pyr']);
set(gca,'xtick',1:length(names),'xticklabel',names,'XTickLabelRotation',20);
set(gca,'box','off','fontsize',15); 

% subplot(2,3,3);
% ok = p(:,8)==1 & abs(d)<1;% & ismember(p(:,6),[2 3]);
% g = [p(ok,1:2),p(ok,6)];% g(:,1) = zscore(g(:,1));
% anovabar(g(:,1:2),g(:,end),'parametric',false);
% ylabel('p(transmission');
% names = {'no fields','field in stim only','field in off only','both'};
% for j=1:4, names{j} = [names{j} ',n=' num2str(sum(g(:,end)==j))]; end
% set(gca,'xtick',1:4,'xticklabel',names,'XTickLabelRotation',20);
% title('pooled CA3 activity -> CA1 int');
% subplot(2,3,3+3);
% g = [abs(d(ok,1)),p(ok,6)]; %g(:,1) = zscore(g(:,1));
% anovabar(g(:,1),g(:,2),'parametric',false);
% ylabel('abs(post-pre)/(post+pre)');
% title('pooled CA3 activity -> CA1 int');
% set(gca,'xtick',1:4,'xticklabel',names,'XTickLabelRotation',20);
% 
% 
% EquateScales(1,2,4);
% EquateScales(3,5,6);