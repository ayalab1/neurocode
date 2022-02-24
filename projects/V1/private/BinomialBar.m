function varargout = BinomialBar(count,n,p)

% Calls the binomial test to plot the binomial test if count/n is different from proportion "p"
% EXAMPLES:
% BinomialBar(sum(heads),sum(heads+tails),0.5); % is the coin biased?
% significant = pValues<0.05; BinomialBar(sum(significant),length(significant),0.05); % are we getting more significant cells than chance levels of p=0.05?
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin == 1 & length(count)==3 % if all 3 inputs were provided in a single vector variable
    p = count(3); n = count(2); count = count(1);
end

ydata = [count/n p]*100;
b = bar(ydata); hold on;

vers = version;
if num2str(vers(1))<8 % Code for Matlab 2010
    bars = get(b,'children'); %where 1st cell is first bars for each row
    for col = 1
        xdata(:,col) = nanmean(get(bars{col},'xdata'));
    end
else
    for col = 1
        xdata(:,col) = b(col).XData + b(col).XOffset;
    end
end
z = zBinomialComparison(count,n,p);

smallnumber = nanstd(ydata(:))/20;
smallnumberx = 0;
p2z = @(p) sqrt(2) * erfcinv(p);;

if abs(z)>1.96, %p<0.025 i.e. 5% two-tailed,
    % draw lines showing which columns are different:
    x1 = xdata(1) + smallnumberx;
    x2 = xdata(2) + smallnumberx;
    [~,this] = max(ydata);
    
    y1 = ydata(this) + smallnumber;
    plot([x1 x2], [y1+smallnumber*sign(y1) y1+smallnumber*sign(y1)], 'k');
    plot([x1 x1], [y1+smallnumber*sign(y1) y1], 'k');
    plot([x2 x2], [y1+smallnumber*sign(y1) y1], 'k');
    % draw stars for the appropriate significance level:
    if abs(z)>p2z(0.0005*2), %p<0.0005 i.e. 0.1% two-tailed
        text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '***','HorizontalAlignment','center');
    elseif abs(z)>p2z(0.005*2), %p<0.005 i.e. 1% two-tailed
        text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '**','HorizontalAlignment','center');
    else
        text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '*','HorizontalAlignment','center');
    end
end

if nargout>0,
    varargout{1} = b;
    varargout{2} = z;
end
