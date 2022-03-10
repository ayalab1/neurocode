function Logical = Unfind(indices, n)
 
%
% When you have a logical vector, 'find' is useful as it gives you the
% indices of the non-zero values.
% This function performs the opposite operation, giving you a logical
% vector with ones in the positions of the indices provided.
%
% If you want the length of the resulting logical to be different from the
% last index provided, provide n (default = last index).
%
 
%% OLD WAY
% 
% Logical = full(sparse(unique(indices), ones(size(unique(indices))), ones(size(unique(indices)))));
% 
% if exist('n', 'var'),
%     if length(Logical)<n,
%         Logical(n) = 0; % add zeroes to the end of the logical so that now it is of length n
%     else
%         Logical = (Logical(1:n)); % if it is shorter, crop it to match n
%     end
% end
% 
% Logical = logical(Logical(:));
% 
%% NEW WAY (17 Feb 2016)
 
if nargin<2,
    n=max(indices);
end
 
Logical = false(n,1);
Logical(indices) = true;