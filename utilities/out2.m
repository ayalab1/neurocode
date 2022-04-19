function output = out2(fun, varargin);

% returns second output argument of fun
% EXAMPLE: order = output2(@sort, matrix, dim);

[~,output] = fun(varargin{:});