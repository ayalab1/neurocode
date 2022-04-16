function varargin = ciplot(x,y,varargin)

%ciplot - alternative to <a href="matlab:help semplot">semplot</a> for medians and confidence intervals
%
% Copyright (C) 2017 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ischar(y),
    varargin = {y,varargin{:}};
    y = x;
    x = 1:size(y,2);
end

color = [0 0 0];
smooth = 0;
portion = [0.05 0.95];

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ciplot">ciplot</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'smooth',
            smooth = varargin{i+1};
            if ~isvector(smooth) || length(smooth) ~= 1,
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help ciplot">ciplot</a>'' for details).');
            end
        case 'color',
            color = varargin{i+1};
%             if ~isvector(color) || length(color) ~= 3,
%                 error('Incorrect value for property ''color'' (type ''help <a href="matlab:help ciplot">ciplot</a>'' for details).');
%                 endcolor = [0 0 0];
%             end
        case 'portion',
            portion = varargin{i+1};
            if any(portion>1),
                warning('You have provided a portion bigger than 1. Percentage is assumed.');
                portion = portion/100;
            end
            if length(portion)==1,
                portion = sort([portion 1-portion]);
            elseif ~isvector(portion) || length(portion) ~= 2,
                error('Incorrect value for property ''portion'' (type ''help <a href="matlab:help ciplot">ciplot</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ciplot">ciplot</a>'' for details).']);
    end
end

if isvector(y),
    handle = plot(x,Smooth(y,smooth),'color',color);
    if nargout>0,
        varargin{1} = handle;
    end
    return
end

xx = [x(:);flipud(x(:))];
yy = [Smooth(quantile(y,portion(1))',smooth); Smooth(flipud(quantile(y,portion(2))'),smooth)];
y = Smooth(nanmedian(y),smooth);

handles = fill(xx,yy,color);
set(handles,'FaceAlpha',0.5,'edgeAlpha',0); 

hold on;
plot(x,y,'color',color,'linewidth',2);
% plot(x,Smooth(nanmedian(y),smooth),'color',color,'linewidth',2);

if nargout>0,
    varargin{1} = handle;
    varargin{2} = a;
end