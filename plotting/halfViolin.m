function varargout = halfViolin(data,varargin)

% data should be arranged in columns; plot half violin for each

% default values
nBins = 1000; 
smooth = 10; 
maxWidth = 0.5; 
colors = {'b','k','g','r','y','m','c'};

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'nbins',
            nBins = varargin{i+1};
        case 'colors'
            colors = varargin{i+1};
        case 'smooth'
            smooth = varargin{i+1}; 
        otherwise,
        error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help anovabox">anovabox</a>'' for details).']);
    end
end

% generate probability distribution
limits = quantile(data(:),[0 1]) + [-1 1]*range(data(:));
[h,ht] = Dist(linspace(limits(1),limits(2),nBins),data);
h = Smooth(h,[smooth 0]);
h = h./max(h(:))*maxWidth;

clear handles
for i=1:size(h,2)
    hh = h(:,i);
    handles(i) = patch( i-1+[hh',zeros(1,numel(ht),1),0],[ht,fliplr(ht),ht(1)],colors{i} ); hold on
    set(handles(i),'EdgeColor','none','FaceAlpha',0.5);
    a = prctile(data(:,i),2.5); b = prctile(data(:,i),97.5);
    plot([i-1,i-1],[a,b],'k','LineWidth',3.0)
    scatter(i-1,nanmean(data(:,i)),800,'k.')
end

end