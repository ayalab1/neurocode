function indices = findmin(data,k)

if nargin<2,
    [~,indices] = min(data);
else
    [~,order] = sort(data,'ascend');
    indices = order(1:k);
end