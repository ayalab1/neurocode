function [indices,values] = findmax(data,k)
%findmax - find local maximum with a  window length of k, if k<2, it calculates global maximum
if nargin<2,
    [values,indices] = max(data);
else
    [values,order] = sort(data,'descend');
    indices = order(1:k);
    values = values(1:k);
end