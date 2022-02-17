function varargout = PETH(data, event, varargin)

% [mat,t,m] = PETH(data,events,<options>).
% 'mat' is a matrix where each row is a peri-event time histogram (PETH) centered on a single event


% default values
duration = [-1 1];
namestring = [inputname(1) ', synchronised to ' inputname(2)];
pictureoptions = {};
smooth = 1;
nBins = 100;
stars = false;
if nargout<1
    show = 'on';
else
    show = 'off';
end
mode = 'l';

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PETH">PETH</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'durations',
            duration = varargin{i+1};
            dt = duration(2);
            if ~isvector(dt) || length(dt) ~= 1,
                error('Incorrect value for property ''dt'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'duration',
            duration = varargin{i+1};
            dt = duration(2);
            if ~isvector(dt) || length(dt) ~= 1,
                error('Incorrect value for property ''dt'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'nbins',
            nBins = varargin{i+1};
            if ~isvector(nBins) || length(nBins) ~= 1,
                error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'show',
            show = varargin{i+1};
            if ~isastring(show,'on','off'),
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'mode',
            mode = varargin{i+1};
            if ~isastring(mode,'l','c'),
                error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'name',
            namestring = varargin{i+1};
            if ~isastring(namestring),
                error('Incorrect value for property ''name'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'smooth',
            smooth = varargin{i+1};
            if ~isvector(smooth) || length(smooth) ~= 1,
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'linetype',
            linetype = varargin{i+1};
            if ~isastring(linetype),
                error('Incorrect value for property ''linetype'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        otherwise,
            pictureoptions{end+1} = varargin{i};%pictureoptions{end+1} = varargin{i+1};
    end
end

seconds_around_event = diff(duration)/2;
% if ~exist('nBins', 'var'),
%     if seconds_around_event < 1,
%         nBins = seconds_around_event*300;
%     else
%         nBins = seconds_around_event*100;
%     end
% end

if size(data,2)==2 % if the provided data is a signal rather than events
    t = linspace(duration(:,1),duration(2),nBins);
    mat_t = bsxfun(@plus,event,t);
    mat = interp1(data(:,1),data(:,2),mat_t);
    m = Smooth(nanmean(mat),smooth);
    if nargout>0,
        varargout{1} = mat;
        varargout{2} = t;
        varargout{3} = m;
    end
    if strcmpi(show,'on'), plot(t', m, pictureoptions{:}); end
    return
end


[sync j] = Sync(data, event, 'durations', duration);
[m, ~, t] = SyncHist(sync, j, 'nBins', nBins, 'smooth', smooth, 'mode', 'mean', 'durations', duration);

if nargout>0,
	s = Bin(sync(:,1),duration,nBins);
	mat = zeros(size(event,1),nBins);
	if isempty(sync);
		t = linspace(duration(1),duration(2),nBins);
		m = zeros(1,nBins);
	else
		if isvector(sync),
			mat(:) = Accumulate(sub2ind(size(mat),j,s),1,numel(mat));
		else
			if strcmp(mode, 'l'); mat(:) = Accumulate(sub2ind(size(mat),j,s),sync(:,end),numel(mat))./Accumulate(sub2ind(size(mat),j,s),1,numel(mat));
			else mat(1:max(sub2ind(size(mat),j,s))) = meanangle(sync(:,end),1,sub2ind(size(mat),j,s)); end
		end
		t = linspace(duration(1),duration(2),nBins);
		% m = mean(mat);
		control = repmat(nanmean(mat(:,~InIntervals(t,7*duration/8)),2),1,nBins);
		[h,p] = ttest(mat,control,0.05/nBins);
	end
	varargout{1} = mat;
	varargout{2} = t;
	varargout{3} = m;
end

if strcmpi(show,'on'),
    if isempty(pictureoptions),
        title([namestring ', ' num2str(numel(j)) ' x ' num2str(numel(unique(j))) ' instances']);
        if exist('linetype','var'), PlotXY(t', m, linetype); else, PlotXY(t', m); end
        if stars, try hold on; plot(t(h>0),m(h>0),'k*'); end; end
        %         try sig = FindInterval(h); PlotIntervals([t(sig(:,1))-mode(diff(t))/2 t(sig(:,2))+mode(diff(t))/2]); end
        title([namestring ', ' num2str(numel(j)) ' x ' num2str(numel(unique(j))) ' instances']);
    else
        title([namestring ', ' num2str(numel(j)) ' x ' num2str(numel(unique(j))) ' instances']);
        PlotXY(t', m, pictureoptions{:}); title([namestring ', ' num2str(numel(j)) ' x ' num2str(numel(unique(j))) ' instances']);
        if stars,try hold on; plot(t(h>0),m(h>0),'k*'); end; end
        %         try sig = FindInterval(h); PlotIntervals([t(sig(:,1))-mode(diff(t))/2 t(sig(:,2))+mode(diff(t))/2]); end
    end
end
