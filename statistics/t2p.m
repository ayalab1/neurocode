function p = t2p(t,df)

%t2p - Transform a t-statistic into a p-value
%
% Input:
%   t   t value
%   df  degrees of freedom
% Output:
%   p   pvalue
%
p = 2*tcdf(-abs(t),df);