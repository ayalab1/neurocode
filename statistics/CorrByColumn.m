function [c,slope,rSquared,SSres] = CorrByColumn(A,B,parametric)

% CorrByColumn - correlate each column of matrix A to each column of matrix B
%
% This is a faster implementation of the loop: 
% for i=1:size(A,1), c(i) = corr(A(:,i),B(:,i)); end
%
% Input: 
%   A           matrix
%   B           matrix
%   parametric  true/false toggle; parametric = true (default) will compute
%               the Pearson coefficient, whereas parametric = 
%               false will compute the Spearman coefficient.
% Output:
%   c           correlation (r)
%   slope       slope of regression
%   rSquared    coefficient of variation (R2)
%   SSres       sum of squares residual
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

bad = isnan(A) | isnan(B);
A(bad) = nan;
B(bad) = nan;
nonconclusive = sum(double(~bad))'<3;

if ~exist('parametric','var'), parametric = true; end

if ~parametric % transform the values to ranks
    A = tiedrank(A);
    B = tiedrank(B);
end

An = bsxfun(@minus,A,nanmean(A,1));
Bn = bsxfun(@minus,B,nanmean(B,1));
An = bsxfun(@times,An,1./sqrt(nansum(An.^2,1)));
Bn = bsxfun(@times,Bn,1./sqrt(nansum(Bn.^2,1)));
c = nansum(An.*Bn,1)';
c(nonconclusive) = nan;

if nargout<2
    return
end
slope = c.*(nanstd(A)./nanstd(B))';
% rSquared = 1-SSres/SStotal; here SStotal = 1
% residual sum of squares is the difference from the fit: sum((An - Bn*c)^2)
rSquared = 1 - sum((An-bsxfun(@times,Bn,c')).^2)';
SSres = (1-rSquared).*(sum(~bad)' - 1).*var(A)';

% Afit = bsxfun(@plus,bsxfun(@times,bsxfun(@minus,B,nanmean(B,1)),slope'),nanmean(A,1));
% startstop = Afit([1 end],:)';
end