function [angle,r] = CircularFun(fun,angles,weights)

%CircularFun - apply any function to circular data
%
%
%  USAGE
%
%    [angle,r] = CircularFun(fun,angles,weights)
%
%    fun            a function handle to the operation to be applied to
%                   the angles
%    angles         vector or matrix of angles (in radians)
%    weights        relative weight given to each angle (default = ones)
%
%  OUTPUT
%
%    angle          resultant angle (in radians)
%    r              resultant vector length
%
%  EXAMPLE
%
% Compute mean phase:
% [m,r] = CircularFun(@mean,phases); 
% This is the same as calling CircularMean
% [s,sr] = CircularFun(@std,phases); 
% would compute the standard deviation of the provided phases
%
% To apply operation along a particular dimention, 
% e.g. to compute the standard deviation along rows of the matrix "phases",
% make a custom inline function:
% fun = @(x) std(x,[],2);
% [s,sr] = CircularFun(fun,phases); 
%
% Copyright (C) 2022 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin<4
    d = find(size(angles)>1,1);
    if isempty(d) %This is a scalar
        m = angles;
        return
    end
end


if ~exist('weights','var') || isempty(weights) || length(weights)==1,
    weights = ones(size(angles));
end

angles = exp(1i*angles).*weights;

complex = 1i*fun(imag(angles)) + fun(real(angles)); % angles in complex form
angle = atan2(imag(complex),real(complex));
r = sqrt(imag(complex).^2+real(complex).^2);

end



