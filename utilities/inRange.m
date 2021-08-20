function out = inRange(x,r,hard)

if hard
    out = x>=r(1) & x<=r(2);
else
    out = x>r(1) & x<r(2);
end