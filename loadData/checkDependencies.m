function [list,matlabVersion] = checkDependencies(filename)

% checkDependencies - Check dependencies for your code

%    =========================================================================
%  USAGE
%
%INPUT
%   [filename]

%    =========================================================================

%OUTPUT
%   [list] - a cell of then full paths to the functions required for your code to run
%   =========================================================================
%
%   EXAMPLE
% list = checkDependencies('getLFP')
% 

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

[list,matlabVersion] = matlab.codetools.requiredFilesAndProducts(filename);
list = list';