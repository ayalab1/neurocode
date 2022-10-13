function [clean, bad, badIntervals] = CleanLFP(lfp,varargin)

% Returns the lfp cleaned from two kinds of artefacts: big artefacts in the
% lfp signal (surpassing a threshold in z-units) and very fast fluctuations,
% in which the signal derivative suprasses a threshold in z-units).
%
%  USAGE
%
%    smoothed = CleanLFP(lfp,<options>)
%
%    data           lfp to clean in [timestamps values] format, where timestamps
%                   are expected to be in secods.
%    <options>      optional list of property-value pairs (see table below))
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  the thresholds for the two signal and derivative artefacts
%                   in z-units (default = [5 10])
%     'aroundArtefact'  time around artefact boundaries that should be removed
%                   (default = [2 0.1]) seconds.
%     'manual'      choose your thresholds manually in debug mode
%    =========================================================================
%
%  OUTPUT
%
%    clean          the lfp in which the rows inside artefacts have been deleted
%    bad            a logical vector indicating which rows of the original lfp file
%                   were flagged as within the artefact and should be removed
%    badIntervals   the [start stop] intervals in which the signal or derivative
%                   surpassed its respective threshold
%
% Dependencies: ConsolidateIntervals, FindInterval
%
% Copyright (C) 2017-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values:
thresholds = [5 10];
aroundArtefact = [0.5 0.1];
manual = false;

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CleanLFP">CleanLFP</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'thresholds',
            thresholds = varargin{i+1};
            if ~isvector(thresholds) || length(thresholds) ~= 2,
                error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help CleanLFP">CleanLFP</a>'' for details).');
            end
        case 'aroundartefact',
            aroundArtefact = varargin{i+1};
            if ~isvector(aroundArtefact) || length(aroundArtefact) ~= 2,
                error('Incorrect value for property ''aroundArtefact'' (type ''help <a href="matlab:help CleanLFP">CleanLFP</a>'' for details).');
            end
        case 'manual',
            manual = varargin{i+1};
            if ~islogical(manual) || length(manual) ~= 1,
                error('Incorrect value for property ''manual'' (type ''help <a href="matlab:help CleanLFP">CleanLFP</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CleanLFP">CleanLFP</a>'' for details).']);
    end
end



% parameters
threshold1 = thresholds(1); % in sigmas deviating from the mean
aroundArtefact1 = aroundArtefact(1); % 2, Big and long artefacts

threshold2 = thresholds(2); % for derivative of z-scored signal
aroundArtefact2 = aroundArtefact(2); % 0.1 Very fast fluctuations (short time scale)

t = lfp(:,1);
values = lfp(:,2);
z = zscore(values);
d = [diff(z);0];
bad = false(size(values));

if manual
    figure; clf; plot(t,z); hold all; plot(xlim,ones(1,2)*threshold1,'r--'); plot(xlim,-ones(1,2)*threshold1,'r--');  ylabel('lfp signal (z-units'); legend('signal','threshold1');
    disp('Feel free to change "threshold1" and type "dbcont" to continue.')
    keyboard;
    if false % Optionally, take a look at the threshold for the derivative
        clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
        disp('Change "threshold2" manually.')
    end
end
% Detect large global artefacts (1)
artefactInterval = t(FindInterval(abs(z)>threshold1));
if numel(artefactInterval)==2
    artefactInterval=artefactInterval(:)';
end
if ~isempty(artefactInterval)
    artefactInterval = ConsolidateIntervals([artefactInterval(:,1)-aroundArtefact1, artefactInterval(:,2)+aroundArtefact1]);
    bad = InIntervals(t,artefactInterval);
else
    artefactInterval = zeros(0,2);
end

%Find noise using the derivative of the zscored signal (2)
noisyInterval = t(FindInterval(abs(d)>threshold2));
if numel(noisyInterval)==2
    noisyInterval=noisyInterval(:)';
end
if ~isempty(noisyInterval)
    noisyInterval = ConsolidateIntervals([noisyInterval(:,1)-aroundArtefact2, noisyInterval(:,2)+aroundArtefact2]);
    bad = bad | InIntervals(t,noisyInterval);
else
    noisyInterval = zeros(0,2);
end

% Substitute noisy signal with interpolated signal as if artefact did not exist
values(bad) = interp1(t(~bad),values(~bad),t(bad,1));

badIntervals = ConsolidateIntervals(sortrows([artefactInterval; noisyInterval]));
ok = ~isnan(values);
clean = [t(ok), values(ok)];

end