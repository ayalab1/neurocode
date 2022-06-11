% batch = StartBatch(@BatchCheeseboardSSA,'OMLcheese.batch');
% X = get(batch,'UserData');
% X15 = X;
% % 
% % %
% % 
% % days = cell2mat(X(:,end-1:end));
% % sesh = repelem((1:size(X,1))',cellfun(@(x) size(x,1),X(:,1)));
% % anglesCell = cat(1,X{:,1}); anglesCellShuffled = cat(1,X{:,2});
% % mm = cell2mat(cellfun(@(x) CircularMean(x(:,1),1,x(:,2))',anglesCell(:,1),'UniformOutput',0));
% % mm0 = cell2mat(cellfun(@(x) CircularMean(x(:,1),1,x(:,2))',anglesCellShuffled(:,1),'UniformOutput',0));
% % ok = ~any(isnan([mm mm0]),2) & days(sesh,2)==1 & days(sesh,1)<3;
% % figure; subplot(2,1,1);
% % anovabar(cos(mm(ok,:))-cos(mm0(ok,:)),days(sesh(ok))+1,'parametric',true); title('15 ms');
% % subplot(2,1,2); anovabar(cos(mm(ok,:))-cos(mm0(ok,:)),days(sesh(ok))+1,'parametric',false); title('medians');
% 
% %%
% batch = StartBatch(@BatchCheeseboardOrderPairs,'OMLcheese.batch');
% X = get(batch,'UserData');
% X50 = X;
% 
% p = cell2mat(X(:,1));
% p0 = cell2mat(X(:,3));
% days = cell2mat(X(:,end-1:end));
% sesh = repelem((1:size(X,1))',cellfun(@(x) size(x,1),X(:,1)));
% 
% figure; 
% d = p-p0; d = d(:,[1 3]); d = diff(d,[],2);
% d = p-p0; d = d(:,[1 2 3]); 
% ok = ~any(isnan(d),2) & days(sesh,2)==1 & days(sesh,1)<3;
% subplot(2,1,1); anovabar(d(ok,:),days(sesh(ok))+1,'parametric',true); title('pyr 50 ms');
% subplot(2,1,2); anovabar(d(ok,:),days(sesh(ok))+1,'parametric',false); title('medians'); 
% 
% % for i=1:8 figure(i); subplot(2,1,1); ylim([-0.01 0.02]); subplot(2,1,2); ylim([-0.01 0.02]); end
% 
%
% %%
% X = cat(1,XsPyrcells{:});
% timescales = [10 15 20 25 50];
% i=3;
% timescale = timescales(i); X = XsPyrcells{i};
% timescale = timescales(i); X = XsAllcells{i};

load('M:\home\raly\temp\X10_X15_X20_X25_X50_BatchCheeseboardOrderPairs.mat','XsAllcells','XsPyrcells');

%%

% X = XppPyr20;
timescales = [10 15 20 25 50];
i=3;
timescale = timescales(i); X = XsPyrcells{i};
y = [-0.05 0.15];
% X = X(21,:);

clf
p = cell2mat(X(:,1));
p0 = cell2mat(X(:,3));
days = cell2mat(X(:,end-1:end));
sesh = repelem((1:size(X,1))',cellfun(@(x) size(x,1),X(:,1)));

% figure; 

d = p-0.5; d = d(:,[1 2 3]); 
ok = ~any(isnan(d),2) & days(sesh,2)==1 & days(sesh,1)<3; d = d*100; y = y*100;
try
subplot(2,2,1); anovabar(d(ok,:),days(sesh(ok))+1,'parametric',true); title(['pyr ' num2str(timescale) ' ms']);
legend('pre','task','post');
legend('location','best','box','off','fontsize',15);
set(gca,'xtick',1:3,'xticklabel',{'sham','MEC','LEC'}); ylim(y);
set(gca,'box','off','tickdir','out','fontsize',15);
ylabel('AB (trials) preference over shuffled data (%)');
end
title(days);
try
subplot(2,2,3); anovabar(d(ok,:),days(sesh(ok))+1,'parametric',false); title('medians'); 
set(gca,'xtick',1:3,'xticklabel',{'sham','MEC','LEC'}); ylim(y);
set(gca,'box','off','tickdir','out','fontsize',15);
ylabel('AB (trials) preference over shuffled data (%)');
legend('pre','task','post');
legend('location','best','box','off','fontsize',15);
end

y = y/2;
d = p; d = d(:,[1 3]); d = diff(d,[],2); d = d*100;
try
    subplot(2,2,2); anovabar(d(ok,:),days(sesh(ok))+1,'parametric',true); title(['pyr ' num2str(timescale) ' ms']);
ylabel('post-pre AB preference (%)');
set(gca,'xtick',1:3,'xticklabel',{'sham','MEC','LEC'}); ylim(y);
set(gca,'box','off','tickdir','out','fontsize',15);
end

subplot(2,2,4); anovabar(d(ok,:),days(sesh(ok))+1,'parametric',false); title('medians'); 
set(gca,'xtick',1:3,'xticklabel',{'sham','MEC','LEC'}); ylim(y);
ylabel('post-pre AB preference (%)');
set(gca,'box','off','tickdir','out','fontsize',15);

































