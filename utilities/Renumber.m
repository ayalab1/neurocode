function x = Renumber(x,mode)

%   This function takes a vector and renumbers all the entries so there are
%   no gaps between them and they start from one.
%   
%   EXAMPLE
%   x = [30 4 4 5 30 6];
%   x = Renumber(x);
%   
%   This returns [1 2 2 3 1 4].


if ~exist('mode','var')
    mode = 'stable';
end

[~,~,x] = unique(x,mode);
