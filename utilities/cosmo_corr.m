function c=cosmo_corr(x,y,corr_type)
% Computes correlation - faster than than matlab's "corr" for Pearson.
%
% c=comso_corr(x[,y[,corr_type]])
%
% Inputs:
%   x          PxM matrix.
%   y          PxN matrix (optional). If omitted then y=x.
%   corr_type  'Pearson' or 'Spearman' or 'Kendall' (optional). If omitted
%              then corrtype='Pearson' and the computation time is
%              significantly reduced for small matrices x and y (with
%              /tiny/ numerical imprecisions) by the use of a custom
%              implementation.
%              Require the matlab stats
% Output:
%   c          MxN matrix with c(i,j)=corr(x(:,i),y(:,j),'type',corr_type).
%
%
% See also: corr
% Copyright: CoSMoMVPA
% aza, sas lab 2017


    y_as_x=false;

    if nargin<2
        corr_type='Pearson';
        y_as_x=true;
    elseif ischar(y)
        corr_type=y;
        y_as_x=true;
    elseif nargin<3
        corr_type='Pearson';
    end

    if y_as_x
        y=x;
    end

    switch corr_type
        case 'Pearson'
            c=corr_pearson(x,y,y_as_x);

        case 'Spearman'
            x_ranks=tiedrank(x);
            y_ranks=tiedrank(y);

            c=corr_pearson(x_ranks,y_ranks,y_as_x);

        otherwise
            % back to Matlab's function (will break if no stat toolbox)
            c=corr(x,y,'type',corr_type);
    end


function c=corr_pearson(x,y,y_as_x)

    % speed-optimized version
    nx=size(x,1);
    ny=size(y,1);

    % subtract mean
    xd=bsxfun(@minus,x,sum(x,1)/nx);
    yd=bsxfun(@minus,y,sum(y,1)/ny);

    % normalization
    n=1/(size(x,1)-1);

    % standard deviation
    xs=(n*sum(xd .^ 2)).^-0.5;
    ys=(n*sum(yd .^ 2)).^-0.5;

    % compute correlations
    c=n * (xd' * yd) .* (xs' * ys);

    if y_as_x
        % ensure diagonal elements are 1
        c=(c+c')*.5;
        dc=diag(c);
        c=(c-diag(dc))+eye(numel(dc));
    end