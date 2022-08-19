function [c,slope,rSquared,SSres,startstop] = CorrByColumn(A,B)

% [c,slope,rSquared] = CorrByColumn(A,B)
% correlate matrices A and B column by column
% i.e. c(i) = corr(A(:,i),B(:,i))
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

Afit = bsxfun(@plus,bsxfun(@times,bsxfun(@minus,B,nanmean(B,1)),slope'),nanmean(A,1));
startstop = Afit([1 end],:)';
end