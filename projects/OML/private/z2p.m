function p = z2p(z)

p = cdf('norm',-(abs(z)),0,1)*2; % 2 as it's two-tailed (we ignore the sign by taking the absolute value)