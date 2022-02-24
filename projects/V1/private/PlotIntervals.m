function varargout = PlotIntervals(intervals, varargin)

% properties:
% facecolor, facealpha, ylim, edgecolor, flip

edge = 'none';
color = [0.7 0.7 0.7];
y = ylim;
lowlim = y(1);
highlim = y(2);
alpha = 0.5;
flip = false;

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case {'facecolor','color'},
            color = varargin{i+1};
        case  {'facealpha','alpha'},
            alpha = varargin{i+1};
        case 'ylim',
            lowlim = varargin{i+1}(1);
            highlim = varargin{i+1}(2);
        case 'edgecolor',
            edge = varargin{i+1};
        case 'flip',
            flip = varargin{i+1};
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).']);
    end
end

if isempty(intervals), return; end




for i=1:size(intervals,1),
    if ~flip,
        h = patch([intervals(i,1) intervals(i,2) intervals(i,2) intervals(i,1)],[lowlim lowlim highlim highlim], color); hold on;
    else
        h = patch([lowlim lowlim highlim highlim],[intervals(i,1) intervals(i,2) intervals(i,2) intervals(i,1)], color); hold on;
    end
    set(h,'edgecolor',edge);
    if alpha~=1,set(h, 'FaceAlpha', alpha);end
    if nargin>0,varargout{1} = h;end
end


if isempty(intervals),
    if ~flip,
        h = patch(zeros(4,1),[lowlim lowlim highlim highlim], color); hold on;
    else
        h = patch([lowlim lowlim highlim highlim],zeros(4,1), color); hold on;
    end
    set(h,'edgecolor','none');
    if alpha~=1,set(h, 'FaceAlpha', alpha);end
    if nargin>0,varargout{1} = h;end
end