function z = zBinomialComparison(s1,n1,s2,n2),

% give count1, n1, count2, n2
% This z-test is equivalent to the chi-square test to see if the probability
% of observing one phenomenon s1/n1, where s1 is the number of positive 
% observations and n1 is the total number of observations, is greater than 
% the probability of observing another phenomenon s2/n2.
% Alternatively, call zBinomialComparison(s1,n1,p) to test if s1/n1 is greater
% than p.
%
% EXAMPLE
% given the logical vector 'significant' 
% zBinomialComparison(sum(significant),length(significant),0.05); 
% tests if the number of significant cases is greater than the number 
% of significant cases to be expected by chance (0.05).
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% if the input was a vector rather then the 4 scalars expected:
if nargin==1 && length(s1)==4,
    n2=s1(4);
    s2=s1(3);
    n1=s1(2);
    s1=s1(1);
end

if nargin<4
%   Only one proportion was given. Will perform Z-test for a proportion (binomial distribution), rather than Z-test for the equality of two proportions
    p = s1/n1;
    if nargin==3, p0 = s2; else p0 = 0.05; end
    n = n1;
    z = (p-p0)./sqrt(p0.*(1-p0)./n);
    return
end

p = [s1./n1 s2./n2];
n = [n1 n2];

P = sum(p.*n,2)./sum(n,2);
z = -diff(p,[],2)./sqrt(P.*(1-P).*(1./n(:,1) + 1./n(:,2)));
