function rho = WeightedCorrCirc(weights,x,alpha)

% Provide a matrix of weights, and this function will check the 
% correlation between the X and Y dimensions of the matrix.
% You can provide the X-values and the (angular) Y-values for 
% each column and row.
%
% Copyright (C) 2017 Céline Drieu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if exist('x','var') && ~isempty('x'),
    if isvector(x), x = repmat(x(:)',size(weights,1),1); end
else
    [x,~] = meshgrid(1:size(weights,2),1:size(weights,1));
end
if exist('alpha','var') && ~isempty('alpha'),
    if isvector(alpha), alpha = repmat(alpha(:),1,size(weights,2)); end
else
    [~,alpha] = meshgrid(1:size(weights,2),1:size(weights,1));
    alpha = 2*pi*alpha/size(weights,1);
end


rxs = WeightedCorr(weights,x,sin(alpha));
rxc = WeightedCorr(weights,x,cos(alpha));
rcs = WeightedCorr(weights,sin(alpha),cos(alpha));

% compute angular-linear correlation (equ. 27.47)
rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));

end

function c = WeightedCorr(weights,x,y)

% Provide a matrix of weights, and this function will check the 
% correlation between the X and Y dimensions of the matrix.
% You can provide the X-values and the Y-values for each column
% and row.
%
% Copyright (C) 2017 Céline Drieu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if exist('x','var') && ~isempty('x'),
    if isvector(x), x = repmat(x(:)',size(weights,1),1); end
else
    [x,~] = meshgrid(1:size(weights,2),1:size(weights,1));
end
if exist('y','var') && ~isempty('y'),
    if isvector(y), y = repmat(y(:),1,size(weights,2)); end
else
    [~,y] = meshgrid(1:size(weights,2),1:size(weights,1));
end

x = x(:);
y = y(:);
w = weights(:);

mX = nansum(w.*x)./nansum(w);
mY = nansum(w.*y)./nansum(w);


covXY = nansum(w.*(x-mX).*(y-mY))./nansum(w(:));
covXX = nansum(w.*(x-mX).*(x-mX))./nansum(w(:));
covYY = nansum(w.*(y-mY).*(y-mY))./nansum(w(:));

c = covXY ./ sqrt(covXX.*covYY);
end