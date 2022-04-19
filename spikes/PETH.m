function varargout = PETH(samples, events, varargin)

%PETH - Compute a peri-event time histogram relative to synchronizing events.
% [mat,t,m] = PETH(data,events,<options>).
% 'mat' is a matrix where each row is a peri-event time histogram (PETH) centered on a single event
%
% Select samples that fall around synchronizing events, and make their
% timestamps relative to the synchronizing events. This can be used to
% build e.g. spike raster plots or successive evoked potentials.
%
%  USAGE
%
%    [matrix,t,mean] = PETH(samples,events,<options>)
%
%    samples        either a list of timestamps (e.g. a spike train) or, in
%                   the case of a continuous signal (e.g. reactivation strength,
%                   local field potential, etc), a matrix of [timestamps values]%
%    events         timestamps to synchronize on (e.g., brain stimulations)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-1 1])
%     'nBins'       number of time bins around the events (default = 101)
%     'mode'        whether the sample data is linear ('l') or circular ('c')
%                   (for example, in the case 'samples' is the phase of an
%                   oscillation).
%     'show'        display the mean PETH (default = 'on' when no outputs are
%                   requested and 'off' otherwise)
%     'smooth'      standard deviation for Gaussian kernel (default = 1 bin)
%                   applied to the mean peri-event activity 'm' (note, no
%                   smoothing will be applied to the output 'matrix')
%     'title'       if the results are displayed ('show' = 'on'), specify a
%                   desired title (default is deduced by variable names)
%     <plot options> any other property (and all the properties that follow)
%                   will be passed down to "plot" (e.g. 'r', 'linewidth', etc)
%                   Because all the following inputs are passed down to "plot",
%                   make sure you put these properties last.
%    =========================================================================
%
%  OUTPUT
%
%    matrix         a matrix containing the counts of a point process (for 
%                   timestamp data) or the avarage activity (for a continous
%                   signal) around the synchronizing events. Each column
%                   corresponds to a particular delay around the event (delay
%                   value indicated in timeBins), and each row corresponds to
%                   a particular instance of "events"
%    timeBins       a vector of time bin delay values corresponding the columns
%                   of 
%    mean           average activity across all events
%
%  EXAMPLE
%
%    PETH(spikes(:,1),stimuli); % show mean spiking activity around the stimuli
%
%    [matrix,timeBins,m] = PETH([lfp.timestamps double(lfp.data(:,1))],deltaWaves.peaks); % compute the mean lfp around delta wave peaks
%    [~,order] = sort(deltaWaves.peakNormedPower); % get the order of delta wave power
%    PlotColorMap(matrix(order,:),'x',timeBins); % plot the mean lfp around delta waves as ordered according to delta wave power
%
%  SEE
%
%    See also Sync, SyncHist, SyncMap, PlotSync, PETHTransition.
%
% Copyright (C) 2018 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% default values
duration = [-1 1];
namestring = [inputname(1) ', synchronised to ' inputname(2)];
pictureoptions = {};
smooth = 1;
nBins = 101;
if nargout<1
    show = 'on';
else
    show = 'off';
end
mode = 'l';

for i = 1:2:length(varargin),
    switch(lower(varargin{i})),
        case 'durations',
            duration = varargin{i+1};
            if ~isvector(duration) || length(duration) ~= 2,
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'duration',
            duration = varargin{i+1};
            if ~isvector(duration) || length(duration) ~= 2,
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
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
        case 'title',
            namestring = varargin{i+1};
            if ~isastring(namestring),
                error('Incorrect value for property ''title'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'smooth',
            smooth = varargin{i+1};
            if ~isvector(smooth) || length(smooth) ~= 1,
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        otherwise,
            pictureoptions = varargin(i:end); break
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

if size(samples,2)==2 % if the provided data is a signal rather than events
    t = linspace(duration(:,1),duration(2),nBins);
    mat_t = bsxfun(@plus,events,t);
    if strcmp(mode, 'l')
        mat = interp1(samples(:,1),samples(:,2),mat_t);
        m = Smooth(nanmean(mat),smooth,'nans','on');
    else % circular data
        unwrapped = unwrap(samples(:,2));
        mat = wrap(interp1(samples(:,1),unwrapped,mat_t));
        angles = nanmean(exp(1i*mat));
        smoothed = [Smooth(imag(angles(:)),smooth,'nans','on') Smooth(real(angles(:)),smooth,'nans','on')];
        m = atan2(smoothed(:,1),smoothed(:,2));
    end
    if strcmpi(show,'on'), plot(t', m, pictureoptions{:}); end
    if nargout>0, varargout{1} = mat; varargout{2} = t; varargout{3} = m; end
    return
else % the samples are a point process
    [sync, j] = Sync(samples, events, 'durations', duration);
    if nargout>0,
        s = Bin(sync(:,1),duration,nBins);
        mat = zeros(size(events,1),nBins);
        if isempty(sync)
            t = linspace(duration(1),duration(2),nBins);
            m = zeros(1,nBins);
        else
            mat(:) = Accumulate(sub2ind(size(mat),j,s),1,numel(mat));
            t = linspace(duration(1),duration(2),nBins);
        end
        varargout{1} = mat;
        varargout{2} = t;
    end
    if strcmpi(show,'on'), % compute 'm' that we will plot
        [m, ~, t] = SyncHist(sync, j, 'nBins', nBins, 'smooth', smooth, 'mode', 'mean', 'durations', duration);
    end
end

if strcmpi(show,'on'),
    if isempty(pictureoptions),
        if exist('linetype','var'), PlotXY(t', m, linetype); else, PlotXY(t', m); end
        title([namestring ', ' num2str(numel(j)) ' x ' num2str(numel(unique(j))) ' instances']);
    else
        PlotXY(t', m, pictureoptions{:}); title([namestring ', ' num2str(numel(j)) ' x ' num2str(numel(unique(j))) ' instances']);
    end
end

if nargout>0,
    varargout{1} = mat;
    varargout{2} = t;
    if nargout>2
        if ~exist('m','var')
            [m, ~, ~] = SyncHist(sync, j, 'nBins', nBins, 'smooth', smooth, 'mode', 'mean', 'durations', duration);
        end
        varargout{3} = m;
    end
end
