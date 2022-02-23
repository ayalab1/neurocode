function z = p2z(p)

% This is two tailed: for p<0.05, it suffices for abs(value)>z
z = sqrt(2) * erfcinv(p);
