function residuals = CorrResiduals(spikes, intervals1, intervals2)

% Computes the residuals of the post vs pre pairwise spike train correlations
% using the intervals provided, as described by Giri2019 JNeurosci
% Pearson's correlation for each pair of neuron's spikes will be computed 
% within intervals1 (e.g. rippes in pre-task sleep) as well as intervals2
% (e.g. ripples in post-task sleep). For each interval, the spikes emitted by each
% unit are counted, and if units M and N fire in a correlated manner within the
% provided windows, corr(M,N) would be high.
%
% SEE CorrInIntervals
%
%
% Copyright (C) 2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

r1 = CorrInIntervals(spikes,intervals1);
r2 = CorrInIntervals(spikes,intervals2);
ok = ~isnan(r1(:)) & ~isnan(r2(:)) & ~tril(ones(size(r1))); % take only cross-correlations (not the diagonal) and take each one once (not t
% linear fit
x = r1(ok); y = r2(ok);
a = polyfit(x(:),y(:),1);
residuals = r2-(r1*a(1)+a(2));

