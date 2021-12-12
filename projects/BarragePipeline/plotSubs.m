function plotSubs(x, y, regID, cellID, UIDs, used, line)
if nargin < 7
    line = '';
end
color = ['#FF00FF';'#77AC30';'#0072BD';'#A2142F';'#4DBEEE';'#000000';'#7E2F8E';'#C233FF'];
type = ['s';'^';'o'];

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(cellID)
    ind = find(UIDs==i,1);
    if ~isempty(ind)
		if ~isempty(find(used==i,1))
			if isempty(line)
				plot(x(ind), y(ind), 'MarkerEdgeColor', color(regID(i),:), 'Marker', type(cellID(i)),'MarkerFaceColor',color(regID(i),:));
			else
				plot(x, y(:,ind), 'MarkerEdgeColor', color(regID(i),:), 'Marker', type(cellID(i)), 'Color', color(regID(i),:),'LineStyle',line,'MarkerFaceColor',color(regID(i),:));
			end
		else
			if isempty(line)
				plot(x(ind), y(ind), 'MarkerEdgeColor', color(regID(i),:), 'Marker', type(cellID(i)));
			else
				plot(x, y(:,ind), 'MarkerEdgeColor', color(regID(i),:), 'Marker', type(cellID(i)), 'Color', color(regID(i),:),'LineStyle',line);
			end
		end
    end
end
h = zeros(11, 1);
h(1) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#77AC30','Color', '#77AC30');
h(2) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#0072BD','Color','#0072BD');
h(3) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#A2142F','Color','#A2142F');
h(4) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#4DBEEE','Color','#4DBEEE');
h(5) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#000000','Color','#000000');
h(6) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#7E2F8E','Color','#7E2F8E');
h(7) = plot(NaN,NaN,'Marker','.','MarkerEdgeColor','#C233FF','Color','#C233FF');
h(8) = plot(NaN,NaN,'k^');
h(9) = plot(NaN,NaN,'ko');
h(10) = plot(NaN,NaN,'ko','MarkerFaceColor','k');
h(11) = plot(NaN,NaN,'Marker','s','MarkerEdgeColor','#FF00FF');
legend(h, 'CA1','CA2','CA3','CTX','DG','MEC','LEC','+Mod','-Mod','BrstDet','Unknown', 'Location', 'eastoutside');
end