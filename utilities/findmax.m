function [indices,values] = findmax(data,k)

if nargin<2,
    [values,indices] = max(data);
else
    [values,order] = sort(data,'descend');
    indices = order(1:k);
    values = values(1:k);
end