function indices = findmin(data,k)
%findmin - find local maximum with a  window length of k, if k<2, it calculates global maximum
if nargin<2,
    [~,indices] = min(data);
else
    [~,order] = sort(data,'ascend');
    indices = order(1:k);
end