function out = gaussian(sigma,numpoints)
%out = gaussian(sigma,numpoints)
%Creates a gaussian curve with area of 1.

x = 1:numpoints;
mu = mean(x);

out = exp( (-(x-mu).^2)./(2*(sigma^2)));
s = sum(out);
out = out/s;

