function [list,matlabVersion] = checkDependencies(filename)

% Check dependencies for your code
%   EXAMPLE
% list = checkDependencies('getLFP')
% "list" is a cell of then full paths to the functions required for your code to run

[list,matlabVersion] = matlab.codetools.requiredFilesAndProducts(filename);
list = list';