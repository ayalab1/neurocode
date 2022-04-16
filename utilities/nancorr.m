function [c,p,n] = nancorr(a,b, varargin)

%nancorr - call <a href="matlab:help corr">corr</a> ignoring nans

% Check if needed:
if nargin>1 && ~any(isnan(a(:))) && ~any(isnan(b(:)))
    [c,p] = corr(a,b,varargin{:});
    if nargout>2, n = size(a,1)*ones(size(c)); end
    return
end
    

if nargin==0 || sum(~isnan(a(:))) == 0 || (nargin>1 && sum(~isnan(b(:))) == 0),
    c = nan;
    p = nan;
    n = nan;
    return
end

try
    if nargin==1,
        c = nan(size(a,2));
        p = nan(size(a,2));
        n = nan(size(a,2));
        for i=1:size(a,2),
            for j=1:size(a,2),
                notnan = ~isnan(a(:,i)) & ~isnan(a(:,j));
                n(i,j) = sum(notnan);
                if n(i,j)>1,
                    [c(i,j) p(i,j)] = corr(a(notnan,i),a(notnan,j));
                end
            end
        end
        
        return
    end
catch
    keyboard;
end

if min(size(a))==1 && min(size(b))==1,
    a = a(:); b = b(:);
    notnan = ~isnan(a(:)) & ~isnan(b(:));
    n = sum(notnan);
    if n<3,
        c=nan;
        if nargout>1,p=nan;end
        return
    end
    if nargout>1,[c,p] = corr(a(notnan), b(notnan),varargin{:});
    else c = corr(a(notnan), b(notnan),varargin{:});
    end
else
    c = nan(size(a,2),size(b,2)); p=c; if nargout>2,n=p;end
    for i=1:size(a,2),
        for j=1:size(b,2),
            [c(i,j),p(i,j)] = nancorr(a(:,i),b(:,j),varargin{:});
            if nargout>2,n(i,j) = sum(~isnan(a(:,i)) & ~isnan(b(:,j))); end
        end
    end
end