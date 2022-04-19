function z = p2z(p)

% transform a z-value (z-units) to a p-value

z = sqrt(2) * erfcinv(p);
