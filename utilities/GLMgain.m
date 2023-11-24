function [gain,shGain,predictions,er,sh,w] = GLMgain(source,target,varargin)

%GLMgain - Train a GLM on your "source" data and see how well it predicts 
% the "target" data. This is measured as a gain improvement relative to
% shuffled versions of your data. To prevent overfitting, the data is 
% divided into a number of sets, and to predict the target data in a given set
% the model is trained on the [source target] data outside of that set.
%
% USAGE
%
%    [gain,gainsShuffledData,predictions,errors,errorsShuffledData,weights] = GLMgain(source,target,<options>)
%
%    source         a vector or a matrix of predictor variables, where each
%                   row is an observation, and each column is a variable
%    target         a vector or a matrix of predicted variables, where each
%                   row is an observation, and each column is a variable. 
%                   The rows of "source" and "target" need to correspond to 
%                   the same observations
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'link'            link function to use in GLM (default = 'log')
%     'dist'            distribution to specify in GLM (default = 'normal',
%                       i.e. minimize devFun = @(mu,y) (y - mu).^2)
%     'nSets'           number of sets the data will be divided into for
%                       the cross-validation step (default = 5)
%     'nIterations'     number of times the algorhithm will estimate the
%                       prediction and the shuffle of the prediction
%     'minEvents'       a minumum number of events that the target variable
%                       is greater than zero in to be predicted (default = 10);
%                       The outputs for variables below this limit will be NaNs
%     'mode'            the way to compute the average error to compute the
%                       gain (either 'mean' or 'median'). Sometimes the errors
%                       are not distributed normally. Using the median guarantees
%                       that shuffled data would produce gains centered on 1,
%                       whereas the mean does not guarantee this (default = 'mean')
%    =========================================================================
%
%   OUTPUT
%
%     gain              estimated gain in prediction quality over the shuffled data
%                       It is equivalent to the average error of the shuffled data 
%                       divided by the average error of the actual data. Values greater 
%                       than 1 indicate better prediction than the shuffled data
%     gainsShuffledData the "gain" values obtained for each of the iterations of the 
%                       shuffle. To get a p-value for your gain, see what proportion
%                       of gainsShuffledData is as great or greater than the gain in 
%                       your data: p = mean(gainsShuffledData>=gain)
%     predictions       a vector or a matrix with dimensions correspoding to the target
%                       variable. predictions(i,j) would therefore be the GLM predicted
%                       value for the j-th target variable based on the values of all
%                       the source variables of the i-th observation.
%     errors            the mean difference between the GLM prediction and the actual data:
%                       errors = average(abs(target - predictions))
%     errorsShuffledData the "errors" values obtained for each of the iterations of the 
%                       shuffle.
%     weights           the coefficients estimated by the GLM. The prediction is obtained
%                       by applying the link function to the source data using these weights.
%                       For example, when using a log-link function: 
%                       prediction = exp(source*weights(2:end) + weights(1))
%                       weights(1) is a constant. Disabling the constant term is not implemented.
%
% Copyright (C) 2020 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values:
nSets = 5;
nIterations = 1000;
minEvents = 10;
link = 'log';
dist = 'poisson';
mode = 'mean';

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help GLMgain">GLMgain</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'link',
            link = varargin{i+1};
        case 'dist',
            dist = varargin{i+1};
        case 'mode',
            mode = varargin{i+1}; 
        case 'nsets',
            nSets = varargin{i+1};
            if ~isscalar(nSets) || mod(nSets,1)>0
                error('Incorrect value for property ''nSets'' (type ''help <a href="matlab:help GLMgain">GLMgain</a>'' for details).');
            end
        case 'niterations',
            nIterations = varargin{i+1};
            if ~isscalar(nIterations) || mod(nIterations,1)>0
                error('Incorrect value for property ''nIterations'' (type ''help <a href="matlab:help GLMgain">GLMgain</a>'' for details).');
            end
        case 'minevents',
            minEvents = varargin{i+1};
            if ~isscalar(minEvents) || mod(minEvents,1)>0
                error('Incorrect value for property ''minEvents'' (type ''help <a href="matlab:help GLMgain">GLMgain</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GLMgain">GLMgain</a>'' for details).']);
    end
end

notempty = sum(target~=0)>minEvents; % at least 10 events
target(:,~notempty) = [];
[nEvents,nUnits] = size(target);
state = warning;
state = state(1).state;
% choose the function (nanmean or nanmedian) that will be used to compute the average error
if strcmp(mode,'mean'), average = @(x) nanmean(x); else average = @(x) nanmedian(x); end
if strcmp(link,'log'), linkFunction = @(data,coefficients) exp(data*coefficients(2:end) + coefficients(1)); 
elseif strcmp(link,'identity'), linkFunction = @(data,coefficients) data*coefficients(2:end) + coefficients(1); 
end
warning('off');

gain = nan(nUnits,1);
shGain = nan(nUnits,nIterations);
er = nan(nUnits,1);
sh = nan(nUnits,nIterations);
w = nan(nUnits,size(source,2)+1,nSets);
predictions = nan(size(target));
for unit = 1:nUnits,
    %% Split data into nSets balanced sets
    setID = nan(nEvents,1);
    zero = target(:,unit)==0;
    ok = false; tic; tictoc = toc;
    while ~ok && ((toc - tictoc) <5) % timeout in 5 seconds to prevent an infinite loop
        % initialise (number sets sequentially)
        setID(zero) = ceil(linspace(1/sum(zero),1,sum(zero))*nSets); % we keep zero and non-zero sets separately for balance
        setID(~zero) = ceil(linspace(1/sum(~zero),1,sum(~zero))*nSets);
        % scramble respective sets sets
        setID(zero) = Scramble(setID(zero));
        setID(~zero) = Scramble(setID(~zero));
        % Make sure there is some variance for each of the sources in each set. Otherwises, re-initialize
        ok = true; for set = 1:nSets, if any(~(max(source(setID~=set,:))>min(source(setID~=set,:)))), ok = false; end; end 
    end
    shuffled = nan(nSets,nIterations);
    errors = nan(nSets,1);
    for set = 1:nSets,
        lastwarn('');
        weights = glmfit(source(setID~=set,:),target(setID~=set,unit),dist,'link',link);
        if strcmp(lastwarn,'Iteration limit reached.') % if the algorhithm failed to converge
            continue
        end
        w(nUnits,:,set) = weights;
        prediction = linkFunction(source(setID==set,:),weights);
        predictions(setID==set,unit) = prediction;
        errors(set,1) = average(abs(target(setID==set,unit) - prediction));
        for iteration = 1:nIterations,
            % get randomly ordered numbers
            [~,indices] = sort(rand(sum(setID==set),1));
            shuffled(set,iteration) = average(abs(target(setID==set,unit) - prediction(indices)));
        end
    end
    er(unit,1) = average(errors(:));
    sh(unit,:) = average(shuffled);
    gain(unit,1) = average(shuffled(:))./average(errors(:));
    shGain(unit,:) = bsxfun(@rdivide,average(shuffled(:)),average(shuffled));
end
warning(state);

if sum(~notempty)>0,
    gain0 = gain;
    shGain0 = shGain;
    er0 = er;
    sh0 = sh;
    w0 = w;
    gain = nan(numel(notempty),1);er = nan(numel(notempty),1);
    sh = nan(numel(notempty),nIterations); shGain = nan(numel(notempty),nIterations);
    w = nan(numel(notempty),size(source,2)+1,nSets);
    gain(notempty) = gain0;er(notempty) = er0;sh(notempty,:) = sh0;
    shGain(notempty,:) = shGain0;
    w(notempty,:,:) = w0;
end


% ------------------------------- Helper function -------------------------------

function x = Scramble(x, n)

% Takes n random samples of x;
%
% Copyright (C) 2018 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if size(x,1) == 1, %if a row vector is provided
    x = x(:); %change it to a column vector
end

if ~exist('n', 'var'),
    n = size(x,1);
end

scr = randperm(size(x,1));
x = x(scr(1:n),:);
