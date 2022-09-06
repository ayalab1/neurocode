%% pool objScore across sessions 
%longer training pool HLR 05/2022
%TO DO add export to csv of pooled results/stats for things importing to py
%TO DO if using matlab to graph change colors for finalization

clearvars;
dirData = 'Y:\OJRproject\';
%dirData = 'Z:\ObjRipple\recordings\';
dirSave = 'Y:\OJRproject\analysis_repo\behavior_longtraining';

animals = {'OJR42','OJR43', 'OJR44', 'OJR45','OJR46','OJR47', 'OJR48','OJR51','OJR52', 'OJR55'}; 


days = {{'day11','day12'};{'day7','day10'};{'day11','day16'};{'day8','day9','day12'};...
        {'day8','day15','day16'};{'day10','day11','day15'};{'day10','day11','day13'};{'day5','day8'};{'day13'};{'day10', 'day12'}};


condition = {[1 1]; [1 1]; [2 3]; [2 3 1]; [2 3 1]; [1 3 2]; [1 3 2]; [2 3]; [3]; [3 2]};


% 1=4h delay 3x training; 2=4h delay 3x training + PFC inh; 3=4h delay 3x training + PFC delay;

%% pool
DIcond = cell(3,1);OPcond = cell(3,1);trainT = cell(3,1);testT = cell(3,1);

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


%% Figure  4h PFC exp
box4hPFC = nan(20,3);
box4hPFC(1:numel(DIcond{1}),1) = DIcond{1};
box4hPFC(1:numel(DIcond{2}),2) = DIcond{2};
box4hPFC(1:numel(DIcond{3}),3) = DIcond{3};

c1=[.23 .24 .23]; %color for plots as RGB - can change for both scatter and boxplot here
c2=[.28 .68 .83]; %color2
c3=[.63 .93 .63];
colors=[c3; c2; c1]; %reverse


figure; 
boxplot(box4hPFC,'Notch','on','Labels',{'control','PFC inh','PFC delay'});hold on;
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
ylabel('Discrimination Index');title('PFC: 3x Training + 4h delayed recall');
set(lines, 'Color', 'k','LineWidth',2);
plot(xlim,[0 0],'--k');hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.4);
    
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
p6=ranksum(DIcond{2},DIcond{3});

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

saveas(gcf,'Y:\OJRproject\analysis_repo\behavior_longtraining\longtraining.fig');
saveas(gcf,'Y:\OJRproject\analysis_repo\behavior_longtraining\longtraining.png');
saveas(gcf,'Y:\OJRproject\analysis_repo\behavior_longtraining\longtraining.pdf');



