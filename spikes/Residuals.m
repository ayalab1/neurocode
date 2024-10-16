function residuals = Residuals(variable, factor)

% Computes the residuals of a variable after removing the part that is
% explained by a factor (second input) using a linear fit.
%
% SEE CorrResiduals
%
%
% Copyright (C) 2024 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

ok = ~isnan(factor(:)) & ~isnan(variable(:));
% linear fit
x = factor(ok); y = variable(ok);
a = polyfit(x(:),y(:),1);
residuals = variable-(factor*a(1)+a(2));

