function z = p2z(p)

%p2z - Transform a p-value into z-units

z = sqrt(2) * erfcinv(p);
