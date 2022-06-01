load('day18_hmmLinearized_2.mat')

%%
l = hmmLinearized.linearized;
t = hmmLinearized.timestamps(:);
edges = hmmLinearized.linearizedBoundaries;

% keep the positions from the first 3 edges (of hmmLinearized.linearizedBoundaries(1:2,:))
% keep the positions from leftbound trials (edges 3:13)
edgesLeft = edges(3:13,:);
edgesRight = edges(14:24,:);

% interpolate positions from rightbound trials as if they were left bound trials:
% fit edges 14:24 to edges 3:13
right = l>edges(14,1);

l(right) = interp1([edgesRight(:,1); edgesRight(end,2)], [edgesLeft(:,1); edgesLeft(end,2)], l(right));

mergedEdges = edges(1:13,:);

%%

angles = l/max(l);