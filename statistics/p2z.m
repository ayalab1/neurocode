function z = p2z(p)

%p2z - Transform a p-value into z-units
%
% Input:
%   p   p value
% Output:
%   z   z unit 
%
z = sqrt(2) * erfcinv(p);
