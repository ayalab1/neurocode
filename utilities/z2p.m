function p = z2p(z)

% Transforms z-value (e.g. the output of zBinomialComparison) into a p-value

p = cdf('norm',-(abs(z)),0,1)*2; % 2 as it's two-tailed (we ignore the sign by taking the absolute value)