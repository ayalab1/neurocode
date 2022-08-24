function p = z2p(z)

%z2p - Transform a z-value into p-value
%
% Input:
%   z   z value
% Output:
%   p   pvalue
%
p = cdf('norm',-(abs(z)),0,1)*2; % 2 as it's two-tailed (we ignore the sign by taking the absolute value)