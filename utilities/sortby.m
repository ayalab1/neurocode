function output = sortby(data,vector2order)

[~,order] = sort(vector2order);
output = data(order,:);