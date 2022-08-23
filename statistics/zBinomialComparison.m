function [z,pValue] = zBinomialComparison(s1,n1,s0,n0)

%zBinomialComparison - Perform a z-test to compare two probabilties.
%
% This z-test is equivalent to the chi-square test to see if the probability
% of observing one phenomenon s1/n1, where s1 is the number of positive
% observations and n1 is the total number of observations, is greater than
% the probability of observing a control phenomenon s0/n0.
% Alternatively, call zBinomialComparison(s1,n1,p) to test if s1/n1 is greater
% than p.
%
%  USAGE
%
%    % To test the difference between two observed probabilities:
%    z = zBinomialComparison(s1,n1,s2,n2)
%
%    s1             number of positive observations for the test phenomenon
%    n1             total number of observations for the test phenomenon
%    s0             number of positive observations for the control phenomenon
%    n0             total number of observations for the control phenomenon
%
%    % To test the difference between one observed probability and a theoretical probability:
%    z = zBinomialComparison(s1,n1,p0)
%
%    s1             number of positive observations for the test phenomenon
%    n1             total number of observations for the test phenomenon
%    p0             theoretical probability of the control phenomenon
%
%  OUTPUT
%
%    z              z-value for the difference between the two probabilities.
%                   For z values exceeding 1.96, the probability of observing
%                   the test phenomenon (s1/n1) is significantly higher than
%                   the probability of observing the control phenomenon. Negative
%                   z-values indicate the test phenomenon is less likely to be
%                   observed than the control phenomenon.
%    pValue         p-value corresponding to the z-value (two-tailed). Note that
%                   two different z-values can result in the same p-value (z=1.96 and
%                   z=-1.96 both result in p=0.05).
%
%  EXAMPLES
%
% % 16 out of 20 interneurons were phase-locked to theta (~80%), compared to only 59 out of 98 pyramidal cells (~60%)
% % Are interneurons significantly more likely to be phase locked than pyramidal cells?
% [z,pValue] = zBinomialComparison(16,20,59,98);  % returns z = 1.6764 and pValue = 0.0937, i.e. the difference is not significant
%
% % An animal was correct in 9 out of 11 trials, is this higher than 50% chance?
% [z,pValue] = zBinomialComparison(9,11,0.5);  % returns z = 2.1106 and pValue = 0.0348 i.e. yes this is significantly more than 50% probability
%
% % I did a statistical test and 26 out of 321 recorded neurons were had a significant p-value at (p<=0.05).
% % That's around 8% of all neurons, is this more than what I should expect by chance (p=0.05)?
% [z,pValue] = zBinomialComparison(26,321,0.05);  % returns z = 2.5481 and pValue = 0.0108 i.e. yes this is significantly more than 5% probability
%
% Copyright (C) 2016-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% if the input was a vector rather then the 4 scalars expected:
if nargin==1 && length(s1)==4
    n0=s1(4);
    s0=s1(3);
    n1=s1(2);
    s1=s1(1);
end

if nargin<4
    %   Only one proportion was given. Will perform Z-test for a proportion (binomial distribution), rather than Z-test for the equality of two proportions
    p = s1/n1;
    if nargin==3
        p0 = s0;
    else
        p0 = 0.05;
    end
    n = n1;
    z = (p-p0)./sqrt(p0.*(1-p0)./n);
    pValue = z2p(z); % transform z-value into a p-value
    return
end

p = [s1./n1 s0./n0];
n = [n1 n0];

P = sum(p.*n,2)./sum(n,2);
z = -diff(p,[],2)./sqrt(P.*(1-P).*(1./n(:,1) + 1./n(:,2)));
pValue = z2p(z); % transform z-value into a p-value
end
% ------------------------------- Helper function -------------------------------

function p = z2p(z)

%z2p - Transform a z-value (e.g. the output of a z-test such as zBinomialComparison) into p-value

p = cdf('norm',-(abs(z)),0,1)*2; % 2 as it's two-tailed (we ignore the sign by taking the absolute value)
end