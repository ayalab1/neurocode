function [clean,isOutlier] = RemoveOutliers(data,whisker)

if ~exist('whisker','var')
    whisker = 1.5;
end
    
% if size(data,2)>0,
%     clean = data; isOutlier = zeros(size(data));
%     for i=1:size(data,2)
%         [clean(:,i),isOutlier(:,i)] = RemoveOutliers(data(:,i));
%     end
% end

q1 = quantile(data,0.25);
q3 = quantile(data,0.75);
d = q3-q1;

isOutlier = bsxfun(@gt,data,q3+whisker*d) | bsxfun(@lt,data,q1-whisker*d);
clean = data;
clean(isOutlier) = nan;