function varargout = semplot(x,y,color,smooth);

% semplot(x,y,color,smooth)

if nargin<2 || (nargin==2 && ischar(y)),
    if nargin==2,
        color = y;
    end
    y = x;
    x = 1:size(y,2);
end

if ~exist('color','var'),
    color = [0 0 0];
end

if ~exist('smooth','var'),
    smooth = 0;
end

if size(y,2)~=length(x),
    y = y';
    if size(y,2)~=length(x),
        error('Y should have one column for each element in X');
    end
end

if isvector(y),
    handles = plot(x,Smooth(y,smooth),'color',color);
    return
end


% a = [nanmean(y)-nansem(y); nanmean(y)+nansem(y)]';
% % a = [quantile(y,0.25); quantile(y,0.75)]';
% a(:,2) = a(:,2)-a(:,1);
% a(:,1) = Smooth(a(:,1),smooth);
% a(:,2) = Smooth(a(:,2),smooth);
% 
% handles = area(x,a,'EdgeColor','none','FaceColor',color);
% delete(handles(1));
% handle = handles(2);
% set(get(handle,'children'),'FaceAlpha',0.5);
% % set(get(handle,'children'),'FaceColor',mean([get(get(handle,'children'),'FaceColor');1 1 1]));


xx = [x(:);flipud(x(:))];
yy = [Smooth(nanmean(y)'-nansem(y)',smooth); Smooth(flipud(nanmean(y)'+nansem(y)'),smooth)];
y = Smooth(nanmean(y),smooth);

handles = fill(xx,yy,color);
set(handles,'FaceAlpha',0.5,'edgeAlpha',0);
% set(handles,'FaceColor',mean([color;1 1 1]),'edgeAlpha',0);

hold on;
plot(x,y,'color',color,'linewidth',2);

if nargout>0,
    varargout{1} = handles;
%     varargout{2} = a;
end