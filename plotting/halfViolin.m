function varargout = halfViolin(data,varargin)

% data should be arranged in columns; plot half violin for each
%
% Copyright (C) 2023 by Can Liu and Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% default values
nColumns = size(data,2);
nBins = 1000; 
binCenters = [];
maxWidth = 0.5; 
colors = get(gca,'ColorOrder'); if size(colors,1)<nColumns, colors = repmat(colors,nColumns,1); end
xdata = 1:nColumns;

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help halfViolin">halfViolin</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'nbins',
            nBins = varargin{i+1};
        case 'bincenters',
            binCenters = varargin{i+1};
        case 'colors'
            colors = varargin{i+1};
        case 'x'
            xdata = varargin{i+1};
        case 'maxwidth'
            maxWidth = varargin{i+1};
        case 'smooth'
            smooth = varargin{i+1};
        case 'grouped'
            data0 = data;
            for j=1:max(data0(:,end)),
                dataCell{j,1} = nan(sum(data0(:,end)==j),max(data0(:,end))); 
                dataCell{j,1}(:,j) = data0(data0(:,end)==j,1:end-1);
            end
            data = cell2mat(dataCell); nColumns = size(data,2); xdata = 1:nColumns;
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help halfViolin">halfViolin</a>'' for details).']);
    end
end
% make sure "colors" is a cell
if isa(colors,'double'), colormatrix = colors; colors = cell(size(colormatrix,1),1); for i=1:size(colormatrix,1), colors{i} = colormatrix(i,:); end; end
if ~isa(colors,'cell'), colors = repmat({colors},nColumns); end
if ~exist('smooth','var'),smooth = nBins/100; end % needs to be defined after "nBins" is provided
if isempty(binCenters), % generate probability distribution
    limits = quantile(data(:),[0 1]) + [-1 1]*range(data(:));
    binCenters = linspace(limits(1),limits(2),nBins);
end

[h,ht] = Dist(binCenters,data);
if nColumns>1,h = Smooth(h,[smooth 0]); else, h = Smooth(h,smooth); end
h = h./max(h(:))*maxWidth;

clear handles
for i=1:nColumns
    hh = h(:,i);
    handles(i) = patch(xdata(i)+[hh',zeros(1,numel(ht),1),0],[ht,fliplr(ht),ht(1)],colors{i} ); hold on
    set(handles(i),'EdgeColor','none','FaceAlpha',0.5);
    a = prctile(data(:,i),2.5); b = prctile(data(:,i),97.5);
    plot(ones(1,2)*xdata(i),[a,b],'k','LineWidth',3.0)
    scatter(xdata(i),nanmean(data(:,i)),800,'k.')
end

if nargout>0,
    varargout{1} = handles;
    varargout{2} = h;
    varargout{3} = ht;
end

