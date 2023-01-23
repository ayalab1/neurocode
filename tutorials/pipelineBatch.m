% Batches are a convenient way to perform pooled analyses

% To compute a batch analysis, first build your batch function, which
% will analyse each session independently.
% This is useful because it makes sure that there is no mix-up between
% sessions. Also, if there is a problem with one session, the analysis
% writes an error message and then moves on to the next session.

% Here, we have already written "pipelineBatch_BatchExample", which you can check
batch = StartBatch(@pipelineBatch_BatchExample,'example.batch');
X = get(batch,'UserData');
% "X" now contains the outputs of the function "pipelineBatch_BatchExample"
% Each row of "X" contains the results of one session. 
% Each column of "X" contains a different output variable of the function "pipelineBatch_BatchExample"
% Here, "pipelineBatch_BatchExample" has three output variables: mPETHs1, mPETHs2, and basepath

 %% pool the results from different sessions:
 % (note that here, all the first outputs of X have the same column number (101) so I can pool them in the same matrix:)
mPETHs1 = cell2mat(X(:,1));
mPETHs2 = cell2mat(X(:,2));
t = linspace(-1,1,101);

subplot(1,3,1);
PlotColorMap(mPETHs1,'x',t);
PlotHVLines(0,'v','k--');
ylabel('Unit ID (pooled)');
xlabel('time from ripple start (s)');
clabel('Hz');
title('pre-task sleep');

subplot(1,3,2);
PlotColorMap(mPETHs2,'x',t);
PlotHVLines(0,'v','k--');
ylabel('Unit ID (pooled)');
xlabel('time from ripple start (s)');
clabel('Hz');
title('post-task sleep');

clims([0 1]*10); % set all clims to [0 10]

subplot(1,3,3);
PlotColorMap(mPETHs2-mPETHs1,'x',t);
PlotHVLines(0,'v','k--');
ylabel('Unit ID (pooled)');
xlabel('time from ripple start (s)');
clabel('Hz difference');
title('difference');
clim([-1 1]*5);


