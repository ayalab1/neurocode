function ppc = PPC(phases,input2)

%PPC - compute Pairwise Phase Consistency across spikes
%
% ppc = PPC(phases) computes the Pairwise Phase Consistency as defined by  
% Vinck, Wingerden, Womelsdorf, Fries, Pennartz (2010, NeuroImage)
% ppc = PPC(phases,trialID) computes P(hat)2 Pairwise Phase Consistency 
% as defined by Vinck, Battaglia, Womelsdorf, Pennartz (2012, J Comp Neurosci)
% ppc = PPC(phases1,phases2) computes the Pairwise Phase Consistency across units
% as defined by Vinck, Wingerden, Womelsdorf, Fries, Pennartz (2010, NeuroImage)
% 
% %  INPUTS
%     phases      a vector containing phases (in radians)
%     input2      it can either be a 
%                 (1) vector of phases (in radians) emitted by another unit
%                 to compute PPC between the two cells, or
%                 (2) a vector of trial IDs (one for each phase). The PPC
%                 will be computed independently for each trial.
% 
%  OUTPUTS
%    ppc      pairwise phase consistency score
%
%
% SEE ALSO
%
% Copyright (C) 2016-2024 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%
%-------------------------------------------------------------------------



if nargin==1
    phases = phases(:);
    difference = bsxfun(@minus,phases,phases');
    difference(eye(size(difference))==1) = nan; % exclude angle comparisons with itself
    % as the diagonal is now composed of nans, there are n*(n-1) non-nan elements in the matrix 'difference'
    % i.e. taking the mean is equilalent to the sum, normalised by n*(n-1) as in Vinck et al (2010)
    ppc = nanmean(cos(difference(:))); % as cos(a-b) = cos(a)*cos(b)+sin(a)*sin(b), Vinck et al's definition of PPC
else
    if isivector(input2) % if it's a vector of integers, assume it's trial ID
        phases = phases(:); trialID = input2(:);
        if numel(trialID)~=numel(phases)
            error('Please provide a trial number for each phase');
        end

        nTrials = max(trialID);
        ppc = nan(nTrials);
        for m=1:nTrials
            for l=1:nTrials
                if l~=m
                    ppc(m,l) = PPC2(phases(trialID==m),phases(trialID==l));
                end
            end
        end
        % as the diagonal is now composed of nans, there are n*(n-1) non-nan elements in the matrix 'ppc'
        % i.e. taking the mean is equilalent to the sum, normalised by n*(n-1) as in Vinck et al (2012)
        ppc = nanmean(ppc(:));
    else % assume it's phases of another unit and compute cross-unit PPC
        phases1 = phases; phases2 = input2;
        difference = bsxfun(@minus,phases1,phases2');
        ppc = nanmean(cos(difference(:))); % as cos(a-b) = cos(a)*cos(b)+sin(a)*sin(b), Vinck et al's definition of PPC
    end
end




