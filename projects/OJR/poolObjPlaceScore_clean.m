%% pool objScore across sessions 
%TO DO add export to csv of pooled results/stats for things importing to py
%TO DO if using matlab to graph change colors for finalization
clearvars;
dirData = 'N:\OJRproject\offline_analysis\';
%dirData = 'Z:\ObjRipple\recordings\';
dirSave = 'N:\OJRproject\analysis_repo\behavior';

animals = {'OJR6','OJR4', 'OJR7', 'OJR5', 'OJR9','OJR10','NG1','NG2','NG3',...
           'NG4','NG5','NG6','OJR11','OJR12','OJR13','OJR18','OJR20','OJR21',...
           'OJR23','OJR24','OJR25','OJR26','OJR27','OJR30','OJR31','OJR32','OJR33','OJR34',...
           'OJR40','MDC7','MGC7'};
days = {{'day2','day4'};{'day3','day9','day10'};{'day1','day4'};{'day3','day4','day7'};{'day2'};...
        {'day1','day3','day4'};{'day9','day7'};{'day8','day9'};
        {'day8','day9'};{'day13','day11'};{'day13','day11'};{'day13','day11'};{'day1','day2','day8'};{'day1'};{'day1'};...
        {'day1','day2','day3','day4'};{'day1','day2','day3','day5','day6'};{'day1','day2','day4','day5','day6'};...
        {'day1','day2'};{'day1','day2','day3','day4'};{'day3','day4','day5'};{'day1','day2'};...
        {'day1','day2','day3','day4','day5','day6'};{'day1','day2','day4','day5','day6','day7'};...
        {'day1','day2','day3','day4','day5'};{'day3','day4'};...
        {'day5','day7','day8','day9','day10','day11'};{'day2','day3','day5','day6','day8'};...
        {'day3','day4'};{'day1027','day1117'};{'day1110','day1118'}};
condition = {[1 2];[2 1 3];[1 2];[2 1 3];[2];[1 3 2];...
            [1 6];[1 6];[1 6];[1 6];[1 6];[1 6];[1 2 3];[1];[1];...
            [4 5 8 7];[4 5 8 6 7];[5 4 8 7 6];...
            [3 1];[1 3 2 6];[7 8 6];[4 6];...
            [7 5 4 8 6 3];[4 5 3 7 8 2];...
            [5 4 7 8 3];[8 7];
            [4 5 7 8 1 6];[5 4 7 8 6];...
            [2 3];[4 5];[4 5]}; 
% 1=4h sham; 2=4h cl; 3=4h ol; 4=4h cl+PFC inh; 5=4h cl+PFC delay;
% 6=1h sham; 7=1h PFC inh; 8=1h PFC delay

%% pool
DIcond = cell(8,1);OPcond = cell(8,1);trainT = cell(8,1);testT = cell(8,1);

for a = 1:length(animals)
    for d = 1:length(days{a})
        cd([dirData animals{a} '\' days{a}{d}]);
        load('objScore.mat');
        
        for c = 1:8
            if condition{a}(d) == c
               DIcond{c} = cat(1,DIcond{c},objScore.discrimination_index);
               %OPcond{c} = cat(1,OPcond{c},objScore.object_preference);               
               trainT{c} = cat(1,trainT{c},objScore.object_training_time(1)+objScore.object_training_time(2));
               testT{c} = cat(1,testT{c},objScore.object_test_time(1)+objScore.object_test_time(2));
            end
        end
        clear objScore;
        
    end
end
%% %%%%%%%%%%%%%%%%% plot       

%% Figure 1: control 
boxCon = nan(20,2);
boxCon(1:numel(DIcond{6}),1) = DIcond{6};
boxCon(1:numel(DIcond{1}),2) = DIcond{1};

c1=[.23 .24 .23]; %color for plots as RGB - can change for both scatter and boxplot here
c2=[.28 .68 .83]; %color2
c3=[.63 .93 .63];

colors=[c2; c1]; %reverse


figure; 
boxplot(boxCon,'Notch','on','Labels',{'1h delay','4h delay'});hold on;
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k','LineWidth',2);
ylabel('Discrimination Index');title('control animals');

plot(xlim,[0 0],'--k');hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.4);
    
end



x=ones(numel(DIcond{6})).*(1+(rand(numel(DIcond{6}))-0.5)/5);
x1=ones(numel(DIcond{1})).*(1+(rand(numel(DIcond{1}))-0.5)/10);
f1=scatter(x(:,1),DIcond{6},'k','filled');f1.MarkerFaceAlpha = 0.8;f1.MarkerFaceColor = c1;f1.MarkerEdgeColor = 'k';hold on 
f2=scatter(x1(:,2).*2,DIcond{1},'k','filled');f2.MarkerFaceAlpha = 0.8;f2.MarkerFaceColor = c2;f2.MarkerEdgeColor = 'k';hold on


p1=signrank(DIcond{6});
p2=signrank(DIcond{1});
p3=ranksum(DIcond{6},DIcond{1});

yt = get(gca,'YTick');  xt = get(gca,'XTick');hold on
axis([xlim floor(min(yt)*1.2) ceil(max(yt)*1.3)])
plot(xt([1 2]), [1 1]*max(yt)*1.15, '-k',  mean(xt([1 2])), max(yt)*1.2);hold on;
text(mean([xt(1),xt(2)]),max(yt)*1.22,['p=' num2str(p3,2)],'FontSize',12);hold on;
text(xt(1),max(yt)*1.05,['p=' num2str(p1,2)],'FontSize',12);hold on;
text(xt(2),max(yt)*1.05,['p=' num2str(p2,2)],'FontSize',12);hold on;
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_control.fig');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_control.png');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_control.pdf');

%% Figure 2: 4h experiment 
box4h = nan(20,3);
box4h(1:numel(DIcond{1}),1) = DIcond{1};
box4h(1:numel(DIcond{2}),2) = DIcond{2};
box4h(1:numel(DIcond{3}),3) = DIcond{3};


colors=[c3; c2; c1]; %reverse

figure; 
boxplot(box4h,'Notch','on','Labels',{'Control','Closed Loop','Delayed'});hold on;
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
ylabel('Discrimination Index');title('4h delayed recall');
set(lines, 'Color', 'k','LineWidth',2);
plot(xlim,[0 0],'--k');hold on;
h2 = findobj(gca,'Tag','Box');
for k=1:length(h2)
    patch(get(h2(k),'XData'),get(h2(k),'YData'),colors(k,:),'FaceAlpha',.4);
    
end



x=ones(numel(DIcond{1})).*(1+(rand(numel(DIcond{1}))-0.5)/5);
x1=ones(numel(DIcond{2})).*(1+(rand(numel(DIcond{2}))-0.5)/10);
x2=ones(numel(DIcond{3})).*(1+(rand(numel(DIcond{3}))-0.5)/15);
f1=scatter(x(:,1),DIcond{1},'filled');f1.MarkerFaceAlpha = 0.8;f1.MarkerFaceColor = c1;f1.MarkerEdgeColor = 'k';hold on 
f2=scatter(x1(:,2).*2,DIcond{2},'filled');f2.MarkerFaceAlpha = 0.8;f2.MarkerFaceColor = c2;f2.MarkerEdgeColor = 'k';hold on
f3=scatter(x2(:,3).*3,DIcond{3},'filled');f3.MarkerFaceAlpha = 0.8;f3.MarkerFaceColor = c3;f3.MarkerEdgeColor = 'k';hold on

p1=signrank(DIcond{1});
p2=signrank(DIcond{2});
p3=signrank(DIcond{3});
p4=ranksum(DIcond{2},DIcond{1});
p5=ranksum(DIcond{3},DIcond{1});

yt = get(gca,'YTick');  xt = get(gca,'XTick');hold on
axis([xlim floor(min(yt)*1.2) ceil(max(yt)*1.4)])
plot(xt([1 2]), [1 1]*max(yt)*1.15, '-k',  mean(xt([1 2])), max(yt)*1.2);hold on;
text(mean([xt(1),xt(2)]),max(yt)*1.22,['p=' num2str(p4,2)],'FontSize',12);hold on;
plot(xt([1 3]), [1 1]*max(yt)*1.25, '-k',  mean(xt([1 3])), max(yt)*1.3);hold on;
text(mean([xt(2),xt(3)]),max(yt)*1.32,['p=' num2str(p5,2)],'FontSize',12);hold on;
text(xt(1),max(yt)*1.05,['p=' num2str(p1,2)],'FontSize',12);hold on;
text(xt(2),max(yt)*1.05,['p=' num2str(p2,2)],'FontSize',12);hold on;
text(xt(3),max(yt)*1.05,['p=' num2str(p3,2)],'FontSize',12);hold on;

saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_4h.fig');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_4h.png');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_4h.pdf');



%% Figure 3: 4h PFC exp
box4hPFC = nan(20,3);
box4hPFC(1:numel(DIcond{1}),1) = DIcond{1};
box4hPFC(1:numel(DIcond{4}),2) = DIcond{4};
box4hPFC(1:numel(DIcond{5}),3) = DIcond{5};




colors=[c3; c2; c1]; %reverse

figure; 
boxplot(box4hPFC,'Notch','on','Labels',{'Control','Closed Loop + PFC Inh','Closed Loop + Delayed'});hold on;
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
ylabel('Discrimination Index');title('PFC: 4h delayed recall');
set(lines, 'Color', 'k','LineWidth',2);
plot(xlim,[0 0],'--k');hold on;
h3 = findobj(gca,'Tag','Box');
for m=1:length(h3)
    patch(get(h3(m),'XData'),get(h3(m),'YData'),colors(m,:),'FaceAlpha',.4);
    
end

x=ones(numel(DIcond{1})).*(1+(rand(numel(DIcond{1}))-0.5)/5);
x1=ones(numel(DIcond{4})).*(1+(rand(numel(DIcond{4}))-0.5)/10);
x2=ones(numel(DIcond{5})).*(1+(rand(numel(DIcond{5}))-0.5)/15);
f1=scatter(x(:,1),DIcond{1},'filled');f1.MarkerFaceAlpha = 0.8;f1.MarkerFaceColor = c1;f1.MarkerEdgeColor = 'k';hold on 
f2=scatter(x1(:,2).*2,DIcond{4},'filled');f2.MarkerFaceAlpha = 0.8;f2.MarkerFaceColor = c2;f2.MarkerEdgeColor = 'k';hold on
f3=scatter(x2(:,3).*3,DIcond{5},'filled');f3.MarkerFaceAlpha = 0.8;f3.MarkerFaceColor = c3;f3.MarkerEdgeColor = 'k';hold on


p1=signrank(DIcond{1});
p2=signrank(DIcond{4});
p3=signrank(DIcond{5});
p4=ranksum(DIcond{4},DIcond{1});
p5=ranksum(DIcond{5},DIcond{1});
p6=ranksum(DIcond{4},DIcond{5});

yt = get(gca,'YTick');  xt = get(gca,'XTick');hold on
axis([xlim floor(min(yt)*1.2) ceil(max(yt)*1.4)])
plot(xt([1 2]), [1 1]*max(yt)*1.15, '-k',  mean(xt([1 2])), max(yt)*1.2);hold on;
text(mean([xt(1),xt(2)]),max(yt)*1.22,['p=' num2str(p4,2)],'FontSize',12);hold on;
plot(xt([1 3]), [1 1]*max(yt)*1.25, '-k',  mean(xt([1 3])), max(yt)*1.3);hold on;
text(mean([xt(2),xt(3)]),max(yt)*1.32,['p=' num2str(p5,2)],'FontSize',12);hold on;
plot(xt([2 3]), [1 1]*max(yt)*1.15, '-k',  mean(xt([2 3])), max(yt)*1.2);hold on;
text(mean([xt(2),xt(3)]),max(yt)*1.22,['p=' num2str(p6,2)],'FontSize',12);hold on;

text(xt(1),max(yt)*1.05,['p=' num2str(p1,2)],'FontSize',12);hold on;
text(xt(2),max(yt)*1.05,['p=' num2str(p2,2)],'FontSize',12);hold on;
text(xt(3),max(yt)*1.05,['p=' num2str(p3,2)],'FontSize',12);hold on;

saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_4h_PFC.fig');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_4h_PFC.png');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_4h_PFC.pdf');



%% Figure 4: 1h PFC exp
box4hPFC = nan(20,3);
box4hPFC(1:numel(DIcond{6}),1) = DIcond{6};
box4hPFC(1:numel(DIcond{7}),2) = DIcond{7};
box4hPFC(1:numel(DIcond{8}),3) = DIcond{8};



colors=[c3; c2; c1]; %reverse


figure; 
boxplot(box4hPFC,'Notch','on','Labels',{'Control','Closed Loop + PFC Inh','Closed Loop + PFC Delayed'});hold on;
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
ylabel('Discrimination Index');title('PFC: 1h delayed recall');
set(lines, 'Color', 'k','LineWidth',2);
plot(xlim,[0 0],'--k');hold on;
h4 = findobj(gca,'Tag','Box');
for n=1:length(h4)
    patch(get(h4(n),'XData'),get(h4(n),'YData'),colors(n,:),'FaceAlpha',.4);
    
end


x=ones(numel(DIcond{6})).*(1+(rand(numel(DIcond{6}))-0.5)/5);
x1=ones(numel(DIcond{7})).*(1+(rand(numel(DIcond{7}))-0.5)/10);
x2=ones(numel(DIcond{8})).*(1+(rand(numel(DIcond{8}))-0.5)/15);
f1=scatter(x(:,1),DIcond{6},'filled');f1.MarkerFaceAlpha = 0.8;f1.MarkerFaceColor = c1;f1.MarkerEdgeColor = 'k';hold on 
f2=scatter(x1(:,2).*2,DIcond{7},'filled');f2.MarkerFaceAlpha = 0.8;f2.MarkerFaceColor = c2;f2.MarkerEdgeColor = 'k';hold on
f3=scatter(x2(:,3).*3,DIcond{8},'filled');f3.MarkerFaceAlpha = 0.8;f3.MarkerFaceColor = c3;f3.MarkerEdgeColor = 'k';hold on



p1=signrank(DIcond{6});
p2=signrank(DIcond{7});
p3=signrank(DIcond{8});
p4=ranksum(DIcond{7},DIcond{6});
p5=ranksum(DIcond{8},DIcond{6});

yt = get(gca,'YTick');  xt = get(gca,'XTick');hold on
axis([xlim floor(min(yt)*1.2) ceil(max(yt)*1.4)])
plot(xt([1 2]), [1 1]*max(yt)*1.15, '-k',  mean(xt([1 2])), max(yt)*1.2);hold on;
text(mean([xt(1),xt(2)]),max(yt)*1.22,['p=' num2str(p4,2)],'FontSize',12);hold on;
plot(xt([1 3]), [1 1]*max(yt)*1.25, '-k',  mean(xt([1 3])), max(yt)*1.3);hold on;
text(mean([xt(2),xt(3)]),max(yt)*1.32,['p=' num2str(p5,2)],'FontSize',12);hold on;
text(xt(1),max(yt)*1.05,['p=' num2str(p1,2)],'FontSize',12);hold on;
text(xt(2),max(yt)*1.05,['p=' num2str(p2,2)],'FontSize',12);hold on;
text(xt(3),max(yt)*1.05,['p=' num2str(p3,2)],'FontSize',12);hold on;

saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_1h_PFC.fig');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_1h_PFC.png');
saveas(gcf,'N:\OJRproject\analysis_repo\behavior\behavior_1h_PFC.pdf');

