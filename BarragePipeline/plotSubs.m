function plotSubs(x, y, regID, cellID, UIDs, used, line)
if nargin < 7
    line = '';
end
color = ['m';'g';'b';'r';'c';'k'];
type = ['s';'^';'o'];

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(cellID)
    ind = find(UIDs==i,1);
    if ~isempty(ind)
		if ~isempty(find(used==i,1))
			if isempty(line)
				plot(x(ind), y(ind), cat(1,color(regID(i)),type(cellID(i))),'MarkerFaceColor',color(regID(i)));
			else
				plot(x, y(:,ind), cat(1,color(regID(i)),type(cellID(i)),line),'MarkerFaceColor',color(regID(i)));
			end
		else
			if isempty(line)
				plot(x(ind), y(ind), cat(1,color(regID(i)),type(cellID(i))));
			else
				plot(x, y(:,ind), cat(1,color(regID(i)),type(cellID(i)),line));
			end
		end
    end
end
h = zeros(8, 1);
h(1) = plot(NaN,NaN,'g');
h(2) = plot(NaN,NaN,'b');
h(3) = plot(NaN,NaN,'r');
h(4) = plot(NaN,NaN,'c');
h(5) = plot(NaN,NaN,'k');
h(6) = plot(NaN,NaN,'k^');
h(7) = plot(NaN,NaN,'ko');
h(8) = plot(NaN,NaN,'ko','MarkerFaceColor','k');
h(9) = plot(NaN,NaN,'sm');
legend(h, 'CA1','CA2','CA3','CTX','DG','+Mod','-Mod','BrstDet','Unknown', 'Location', 'eastoutside');
end