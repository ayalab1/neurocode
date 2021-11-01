function plotSubs(x, y, regID, cellID, UIDs, line)
if nargin < 6
    line = '';
end
color = ['m';'g';'b';'r';'c';'k'];
type = ['*';'^';'o'];

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(cellID)
    ind = find(UIDs==i,1);
    if ~isempty(ind)
        if isempty(line)
            plot(x(ind), y(ind), cat(1,color(regID(i)),type(cellID(i))));
        else
            plot(x, y(:,ind), cat(1,color(regID(i)),type(cellID(i)),line));
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
h(8) = plot(NaN,NaN,'*m');
legend(h, 'CA1','CA2','CA3','CTX','DG','+Mod','-Mod','Unknown', 'Location', 'eastoutside');
end