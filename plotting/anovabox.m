function varargout = anovabox(data,groups,varargin)

%anovabox - Plots box-and-whisker plot for each group
%
%   Plots box-and-whisker plot for each group, indicating significant
%   differences between the groups. Each group is also tested against 0.
%   For non-parametric data, median +/- standard error of the median are 
%   displayed instead, and non-parametric tests are performed.
%
%  USAGE
%
%    [handles,h] = anovabox(data,groups,<options>)
%
%    data           data to test, provided either in matrix (each column is
%                   a different group) or a vector (grouping indicated through
%                   a grouping variable; see examples below)
%    group          grouping variable, indicating the group that each of the
%                   rows of 'data' corresponds to. Alternatively, it can be
%                   empty (default; only possible if 'data' is a matrix) or
%                   a string 'group' (indicating that the grouping variable
%                   is the last column of 'data')
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%    'alpha'        two values for the confidence levels for the following 
%                   tests: (1) between the data and 0, and (2) between
%                   groups (default = [0.05 0.05]). Set alpha to zero to 
%                   skip respective test (e.g. alpha = [0 0.05] would skip
%                   tests between each column and zero).
%    'parametric'   either 'on' (to display means and perform t-tests) 
%                   or 'off' (to display median and perform non-parametric 
%                   Wilcoxon signed-rank and ranksum tests / friedman).
%    'correction'   type of statistical correction to use to control for
%                   multiple comparisons, i.e. between all pairs of multiple
%                   groups (default = tukey-kramer)
%    'precedence'   in case of a two-way anova, indicate which grouping variable
%                   should be tested while the other variable remains stable
%                   By default (default = 1), the first variable, i.e. the columns 
%                   of the matrix 'data' will be tested for a difference in 
%                   each of the groups described by the 2nd grouping variable.
%                   Set to 2 to inverse this behavior.                 
%    'paired'       a setting indicating whether values in the same row are
%                   paired (default = 'on').
%    'lock'         prevent the function from changing the y-axis limits
%                   (default = 'off').
%    'boxwidth'     the width of the boxes to be plotted (default = 0.3).
%    'ns'           display the "n.s." label when between-group comparisons
%                   are not significant (default = false).
%    =========================================================================
%
%  OUTPUT
%
%    handles        a structure of handles for the objects in the figure.
%                   Includes the fields 'bar','errorbars','stars' (for the 
%                   stars above each bar indicating that values in that group
%                   are significantly different from 0), and 'comparisons'
%                   for all the objects (lines and stars) illustrating the 
%                   between-group comparisons.
%    h              h > 0 if the null hypothesis (that the groups are not 
%                   different) can be rejected (h values of 1, 2, and 3, 
%                   correspond to confidence levels 0.05, 0.01, and 0.001, 
%                   respectively)
%
%   Provide data as a matrix (each column will be treated separately), or
%   as grouped data with "group" as a grouping vector (each line indicating
%   which group the respective row in "data" corresponds to). Alternatively,
%   the grouping vector could be the last column of 'data' itself, in which case
%   you can call anovabar(data,'grouped');
%
%  EXAMPLES
%
%  Case 1: the two columns of 'data' correspond to two paired conditions
%  (e.g. control firing 2s before the intervals of interest and
%  the response of interest). Each row corresponds to paired observations:
%       data = [CountInIntervals(spikes(:,1),intervals-2) CountInIntervals(spikes(:,1),intervals)];
%       anovabar(data); % or anovabar(data,[]);
%  Note that for parametric data, the anova will not take the pairing into 
%  account (but for non-parametric data, a paired friedman test is performed).
%
% Case 2: the observations are not paired; groups are therefore indicated
% separately by a grouping variable:
%       controlData = CountInIntervals(spikes(:,1),baseline);
%       % note that here, 'baseline' intervals are not necessarily paired to the intervals of interest below
%       responseData = CountInIntervals(spikes(:,1),intervals);
%       data = [controlData; responseData]; groups = [ones(size(controlData)); ones(size(responseData))*2];
%       anovabar(data,groups);
%
% Case 2bis: same as Case 2, but the grouping variable is contained within 'data':
%       controlData = CountInIntervals(spikes(:,1),baseline); controlData(:,end+1) = 1;
%       responseData = CountInIntervals(spikes(:,1),intervals); responseData(:,end+1) = 2;
%       data = [controlData; responseData];
%       anovabar(data,'grouped'); % or, of course, anovabar(data(:,1),data(:,end)) as in case 2
%
% Case 3: two-way anova: data are grouped according to the columns in 'data'
% AND ALSO by the grouping variable. In this case, the difference between the
% (paired) columns will be tested for each of the groups indicated by the grouping
% variable. For the opposite behavior (testing for the difference between the
% (unpaired) groups indicated by the grouping variable, for each of the columns
% of 'data', use <a href="matlab:help anovabar2">anovabar2</a>.
%
% Copyright (C) 2018-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
alpha = [0.05 0.05];
correction = 'tukey-kramer'; 
grouped = false;
data = double(data);
parametric = true;
precedence = 1;
paired = true;
lock = false; % lock y axis
boxwidth = 0.3;
ns = false;

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'alpha',
            alpha = varargin{i+1};
            if ~isdvector(alpha) || length(alpha)>2
                error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'parametric',
            parametric = varargin{i+1};
            if isastring(lower(parametric),'on','off')
                if strcmpi(parametric,'on'); parametric = true; else, parametric = false; end
            end
            if ~islogical(parametric) || length(parametric)>1
                error('Incorrect value for property ''parametric'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'paired',
            paired = varargin{i+1};
            if isastring(lower(paired),'on','off')
                if strcmpi(paired,'on'); paired = true; else, paired = false; end
            end
            if ~islogical(paired) || length(paired)>1
                error('Incorrect value for property ''paired'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'correction'
            correction = varargin{i+1};
            if ~isastring(correction)
                error('Incorrect value for property ''correction'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'precedence'
            precedence = varargin{i+1};
            if ~isdvector(precedence,'#1'),
                error('Incorrect value for property ''precedence'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'boxwidth'
            boxwidth = varargin{i+1};
            if ~isdvector(boxwidth,'#1'),
                error('Incorrect value for property ''boxwidth'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'ns'
            ns = varargin{i+1};
            if isastring(lower(ns),'on','off')
                if strcmpi(ns,'on'); ns = true; else, ns = false; end
            end
            if ~islogical(ns) || length(ns)>1
                error('Incorrect value for property ''ns'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        case 'lock'
            lock = varargin{i+1};
            if isastring(lower(lock),'on','off')
                if strcmpi(lock,'on'); lock = true; else, lock = false; end
            end
            if ~islogical(lock) || length(lock)>1
                error('Incorrect value for property ''lock'' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).']);
    end
end

% if a single value for alpha is provided, use that value for both kinds of comparisons (to zero and between groups):
if length(alpha)==1, alpha = [alpha alpha]; end
if parametric
    test0 = @(x) out2(@ttest,x);
    testbetween =  @anova1;
    if paired, testpaired = @(x) anova2(x,1,'off'); else testpaired = @(x) testbetween(x,[],'off'); end
else
    test0 = @(x) signrank(x);
    testbetween = @kruskalwallis;
    if paired,
        if size(data,2)==2, testpaired = @helper_signrank; else, testpaired = @(x) friedman(x,1,'off'); end
    else testpaired = @(x) testbetween(x,[],'off');
    end
end
if nargin==1 || isempty(groups)
    groups = [];
    if sum(sum(~isnan(data))==0) > 0, %if one of the groups contains no real values
        return
    end
else
    grouped = true;
    if isastring(groups,'groups','grouped','group')
        groups = data(:,end);
        data = data(:,1:end-1);
    end
    ok = ~isnan(groups);
    data = data(ok,:);
    groups = groups(ok);
    groups = double(groups);
end

nGroups = numel(unique(groups));
if ~isempty(groups)
    [g1,g2] = meshgrid(1:size(data,2),1:size(data,1));
    if precedence==1 % separate the groups provided by the grouping variable; test difference between columns of 'data' in the same group
        hbox = boxchart(groups(g2(:)),data(:),'MarkerStyle','none','GroupByColor',g1(:),'boxwidth',boxwidth);
        u = unique(groups);
        xData = bsxfun(@plus,u,linspace(-0.5 + (1/(size(data,2)*2)),0.5 - (1/(size(data,2)*2)),size(data,2)));
        for j=1:size(data,2), for i=1:nGroups
                clean = RemoveOutliers(hbox(j).YData(hbox(j).XData==u(i)));
                try yData01(i,j) = min(clean); yData02(i,j) = max(clean);
                end; end; end
        yData2 = yData02 + 0.2*sqrt(std([yData01(:);yData02(:)]));
        yData1 = yData01 - 0.2*sqrt(std([yData01(:);yData02(:)]));
        yData = nanmean(cat(3,yData01,yData02),3);
        if ~lock, try ylim([min([yData1(:);yData2(:)]) max([yData1(:);yData2(:)])]); end; end
        yData2 = yData02 + 0.1*sqrt(std([yData01(:);yData02(:)]));
    else
        % separate the columns of 'data'; test difference between groups provided by the grouping variable in the same column
        hbox = boxchart(g1(:),data(:),'MarkerStyle','none','GroupByColor',groups(g2(:)),'boxwidth',boxwidth);
        xData = g1(1:nGroups,:);
        xData = bsxfun(@plus,xData,linspace(-0.5 + (1/(nGroups*2)),0.5 - (1/(nGroups*2)),nGroups)');
        for j=1:nGroups, for i=1:size(data,2),
                clean = RemoveOutliers(hbox(j).YData(hbox(j).XData==i));
                try yData01(j,i) = min(clean); yData02(j,i) = max(clean);
                end; end; end
        yData2 = yData02 + 2*sqrt(std([yData01(:);yData02(:)]));
        yData1 = yData01 - 2*sqrt(std([yData01(:);yData02(:)]));
        yData = nanmean(cat(3,yData01,yData02),3);
        if ~lock, ylim(mmax([yData1(:);yData2(:)])); end
        yData2 = yData02 + sqrt(std([yData01(:);yData02(:)]));
    end
else
    [g1,g2] = meshgrid(1:size(data,2),1:size(data,1));
    hbox = boxchart(g1(:),data(:),'MarkerStyle','none','boxwidth',boxwidth);
    % add x and yData!
     xData = 1:size(data,2);
     for i=1:size(data,2), clean = RemoveOutliers(hbox.YData(hbox.XData==i));   try yData01(i) = min(clean); yData02(i) = max(clean); end; end
     yData2 = yData02 + 2*sqrt(std([yData01(:);yData02(:)]));
     yData1 = yData01 - 2*sqrt(std([yData01(:);yData02(:)]));
     yData = nanmean(cat(3,yData01,yData02),3);
     if ~lock, try ylim([min([yData1(:);yData2(:)]) max([yData1(:);yData2(:)])]); end; end
     yData2 = yData02 + sqrt(std([yData01(:);yData02(:)]));
end
hold on;

%% Adding significance indication (stars etc)
nStars0 = 0; nObjectsBetween = 0;
% one way anova
if ~grouped || size(data,2)==1 
    if grouped, [~,~,stats] = testbetween(data(:,1),groups,'off'); u = unique(groups)';
    else [~,pTest,stats] = testpaired(data); u = 1:size(data,2);
    end
    if alpha(1)>0 % Unless alpha(1)=0, for each group, check if it's different from zero
        for i=1:length(u)
            if grouped, thesedata = data(groups==u(i),1); else thesedata = data(:,i); end
            if sum(~isnan(thesedata))>2
                p = test0(thesedata);
                if p<alpha(1)
                    % Put a little star above it to show it's significantly different from zero
                    nStars0 = nStars0+1;
                    if sign(yData(i))>0, thisy = mean([max(yData2(i)),max(yData02(i))]);
                    else, thisy = mean([min(yData1(i)),min(yData01(i))]); 
                    end
                    handlesStars0(nStars0) = plot(xData(i), thisy, 'k*', 'markersize', 5, 'MarkerEdgeColor', [0 0 0]);
                end
            end
        end
    end
    u = unique(xData)';
    sigFor2Groups = nan(length(u),1);
    if alpha(2)>0 % Unless alpha(2)=0, perform tests between groups
        if ischar(stats) && strcmp(stats,'signrank') % signrank test doesn't give a "stats" output compatible with multcompare
            comparison = [1 2 (pTest<alpha(2)) + (pTest<alpha(2)/5) (pTest<alpha(2)/50)];
        else
            comparison = multcompare(stats,'display', 'off','alpha',alpha(2),'ctype',correction); % default alpha = 0.05
            comparison = [comparison(:,1:2) double(comparison(:,3).*comparison(:,5)>0)]; %if the upper and lower bound have the same sign
            comparison2 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/5,'ctype',correction); % **, default alpha = 0.01
            comparison(:,3) = comparison(:,3) + double(comparison2(:,3).*comparison2(:,5)>0); % third column shows the number of stars to be included. 1 for 0.05, 1 more for 0.01, and another one for 0.001       if sum(double(comparison2(:,3).*comparison2(:,5)>0)),
            comparison3 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/50,'ctype',correction); % ***, default alpha = 0.001
            comparison(:,3) = comparison(:,3) + double(comparison3(:,3).*comparison3(:,5)>0);
        end
        % Now to plot them:
        for i=1:size(comparison,1)
            s = sign(Portion(yData(:)>0)./(Portion(yData(:)>0 | yData(:)<0))-0.5); if s==0,s=sign(nanmean(yData(:)));end
            smallnumber = 0.1*diff(ylim); %for display purposes, so things don't overlap
            y1 = s*max(yData2(comparison(i,1:2)))+smallnumber/2*((comparison(i,2)-comparison(i,1)-0.5));
            y2 = y1+smallnumber*0.1;
            % make sure the stars are not outside the axis
            if y2>max(ylim), y2 = max(ylim); y1 = y2-smallnumber*0.05; end
            x1 = xData(comparison(i,1));
            x2 = xData(comparison(i,2));
            if comparison(i,3)==0
                if ns, text(mean([x1 x2]), y2, 'n.s.'); end
            else
                plot([x1 x2], [y1 y1], 'k');
                plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                text(mean([x1 x2]), y2, repmat('*',1,comparison(i,3)));
            end
        end
    end
    comparison = sigFor2Groups;

% two way anova
else 
    u = unique(groups)';
    % Test if each group is different from zero and draw stars
    if alpha(1)>0 % If alpha(1)=0, skip this
        for j=1:length(u), for i=1:size(data,2)
                if precedence==1, thisx = xData(j,i); thisy = sign(yData(j,i))*mean([yData2(j,i) yData02(j,i)]);
                else thisx = xData(j,i); thisy = sign(yData(j,i))*(mean([yData2(j,i) yData02(j,i)]) + nanmean(abs(yData(:)))/10);
                end
                thesedata = data(groups==u(j),i);
                if sum(~isnan(thesedata))>2, p = test0(thesedata);
                    if p<alpha(1), nStars0 = nStars0+1;
                        handlesStars0(nStars0) = plot(thisx, thisy, 'k*', 'markersize', 5, 'MarkerEdgeColor', [0 0 0]);
                    end % Put a little star above it
                end
            end; end
    end
    % Test if groups are different between each other
    sigFor2Groups = nan(length(u),1);
    if alpha(2)>0 % If alpha(2)=0, skip this
        if precedence==1, maxj = length(u); else, maxj = size(data,2); end
        for j=1:maxj
            if precedence==1, [~,~,stats] = testpaired(data(groups==u(j),:));
            else, [~,~,stats] = testbetween(data(:,j), groups,'off');
            end
            try
        		comparison = multcompare(stats,'display', 'off','alpha',alpha(2),'ctype',correction);
            catch
                continue % if a whole column of data is missing, then comparisons are impossible
            end
            comparison = [comparison(:,1:2) double(comparison(:,3).*comparison(:,5)>0)]; %if the upper and lower bound have the same sign
            comparison2 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/5,'ctype',correction);
            comparison(:,3) = comparison(:,3) + double(comparison2(:,3).*comparison2(:,5)>0); % third column shows the number of stars to be included. 1 for 0.05, 1 more for 0.01, and another one for 0.001       if sum(double(comparison2(:,3).*comparison2(:,5)>0)),
            comparison3 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/50,'ctype',correction);
            comparison(:,3) = comparison(:,3) + double(comparison3(:,3).*comparison3(:,5)>0);
            sigFor2Groups(j,1) = comparison(1,3);
            for i=1:size(comparison,1)
                smallnumber = 0.1*diff(ylim); %for display purposes, so things don't overlap
                if precedence==1, 
                    s = sign(nanmean(nanmean(data(groups==u(j),:)))); 
                    y1 = s*max(yData2(j,:))+smallnumber/3*(comparison(i,2)-comparison(i,1));
                    x1 = xData(j,comparison(i,1));
                    x2 = xData(j,comparison(i,2));
                else, 
                    s = sign(nanmean(data(:,j))); 
                    y1 = s*max(yData2(:,j))+smallnumber/3*(comparison(i,2)-comparison(i,1));
                    x1 = xData(comparison(i,1),j);
                    x2 = xData(comparison(i,2),j);
                end
                y2 = min([y1+smallnumber,max(ylim)]);
                if x1>x2, this = x1; x1 = x2; x2 = x1; end
                smallnumberx = diff(sort(xData(:)))/10; smallnumberx = smallnumberx(1);
                x1 = x1 + smallnumberx; x2 = x2 - smallnumberx;
                if comparison(i,3)==1
                    plot([x1 x2], [y1 y1], 'k');
                    plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                    plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                    text(mean([x1 x2]), y2, '*','HorizontalAlignment','center'); 
                elseif comparison(i,3)==2
                    plot([x1 x2], [y1 y1], 'k');
                    plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                    plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                    text(mean([x1 x2]), y2, '**','HorizontalAlignment','center');
                elseif comparison(i,3)==3
                    plot([x1 x2], [y1 y1], 'k');
                    plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                    plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                    text(mean([x1 x2]), y2, '***','HorizontalAlignment','center');
                end
            end
        end
    end
    comparison = sigFor2Groups;
end

hold off

if nargout>0
    if nStars0==0, handlesStars0 = []; end
    if nObjectsBetween==0, handlesComparisons = []; end
    handles.bar = hbox;
    handles.errorbars = hErrorbar;
    handles.stars = handlesStars0;
    handles.comparisons = handlesComparisons;
    varargout{1} = handles;

    varargout{2} = comparison;
end

% ------------------------------- Helper functions -------------------------------

function [h,p,stats] = helper_signrank(data,varargin)

ranks = data; ranks(:) = tiedrank(data(:));
p = signrank(diff(data,[],2));
h = p<0.05;
stats = 'signrank';



