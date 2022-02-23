function portion = Portion(logical)

if isvector(logical), logical = logical(:); end

portion = sum(logical)/size(logical,1);