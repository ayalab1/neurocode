function [W, H, cost,loadings,power] = seqNMF(X, varargin)

%seqNMF - Factorizes the NxT data matrix X into K factors
% USAGE: 
%
% [W, H, cost, loadings, power] = seqNMF(X, ...    % X is the data matrix
%       'K', 10, 'L', 20, 'lambda', .1, ...        % Other inputs optional
%       'W_init', W_init, 'H_init', H_init, ...
%       'showPlot', 1, 'maxiter', 20, 'tolerance', -Inf, 'shift', 1, ... 
%       'lambdaL1W', 0, 'lambdaL1H', 0, ...
%       'lambdaOrthoH', 0, 'lambdaOrthoW', 0, 'M', M)
% DESCRIPTION:
%
%   Factorizes the NxT data matrix X into K factors 
%   Factor exemplars are returned in the NxKxL tensor W
%   Factor timecourses are returned in the KxT matrix H
%
% CREDITS: Emily Mackevicius and Andrew Bahle, 2/1/2018
% https://www.biorxiv.org/content/early/2018/03/02/273128
% Adapted to neurocode by aza, 2021
%
%
% INPUTS:
%
%  X                                                    Data matrix (NxT) to factorize
% 'K'               10                                  Number of factors
% 'L'               100                                 Length (timebins) of each factor exemplar
% 'lambda'          .001                                Regularization parameter
% 'W_init'          max(X(:))*rand(N,K,L)               Initial W
% 'H_init'          max(X(:))*rand(K,T)./(sqrt(T/3))    Initial H (rows have norm ~1 if max(data) is 1)
% 'showPlot'        1                                   Plot every iteration? no=0
% 'maxiter'         100                                 Maximum # iterations to run
% 'tolerance'       -Inf                                Stop if improved less than this;  Set to -Inf to always run maxiter
% 'shift'           1                                   Shift factors to center; Helps avoid local minima
% 'lambdaL1W'       0                                   L1 sparsity parameter; Increase to make W's more sparse
% 'lambdaL1H'       0                                   L1 sparsity parameter; Increase to make H's more sparse
% 'W_fixed'         0                                   Fix W during the fitting proceedure   
% 'SortFactors'     1                                   Sort factors by loadings
% 'lambdaOrthoH'    0                                   ||HSH^T||_1,i~=j; Encourages events-based factorizations
% 'lambdaOrthoW'    0                                   ||Wflat^TWflat||_1,i~=j; ; Encourages parts-based factorizations
% 'useWupdate'      1                                   Wupdate for cross orthogonality often doesn't change results much, and can be slow, so option to remove  
% 'M'               ones(N,T)                           Masking matrix if excluding a random test set from the fit
% ------------------------------------------------------------------------
% OUTPUTS:
%
% W                         NxKxL tensor containing factor exemplars
% H                         KxT matrix containing factor timecourses
% cost                      1x(#Iterations+1) vector containing 
%                               reconstruction error at each iteration. 
%                               cost(1) is error before 1st iteration.
% loadings                  1xK vector containing loading of each factor 
%                               (Fraction power in data explained by each factor)
% power                     Fraction power in data explained 
%                               by whole reconstruction
%
%                           Note, if doing fit with masked (held-out) data,
%                               the cost and power do not include masked
%                               (M==0) test set elements
% ------------------------------------------------------------------------

% Parse inputs 
p=inputParser;
addParameter(p,'spikes',{},@isstruct);
addParameter(p,'Interval',[],@isnumeric);
addParameter(p,'bin',20e-3,@isnumeric);
addParameter(p,'window',20e-3,@isnumeric);

addParameter(p,'K',10);
addParameter(p,'L',100);
addParameter(p,'lambda',.001);
addParameter(p,'showPlot',1);
addParameter(p,'maxiter',100);
addParameter(p,'tolerance',-Inf);
addParameter(p,'shift',1);
addParameter(p,'lambdaL1W',0);
addParameter(p,'lambdaL1H',0);
addParameter(p,'W_fixed',0);
addParameter(p,'W_init', nan); % depends on K--initialize post parse
addParameter(p,'H_init', nan); % depends on K--initialize post parse
addParameter(p,'SortFactors', 1); % sort factors by loading?
addParameter(p,'lambdaOrthoW',0); % for this regularization: ||Wflat^TWflat||_1,i~=j
addParameter(p,'lambdaOrthoH',0); % for this regularization: ||HSH^T||_1,i~=j
addParameter(p,'useWupdate',1); % W update for cross orthogonality often doesn't change results much, and can be slow, so option to skip it 
addParameter(p,'M',nan); % Masking matrix: default is ones; set elements to zero to hold out as masked test set


parse(p,varargin{:});
spikes = p.Results.spikes;
Interval = p.Results.Interval;
bin = p.Results.bin;
window = p.Results.window;

K = p.Results.K; 
L = p.Results.L; 
lambda = p.Results.lambda; 
showPlot = p.Results.showPlot; 
maxiter = p.Results.maxiter; 
tolerance = p.Results.tolerance; 
shift = p.Results.shift; 
lambdaL1W = p.Results.lambdaL1W; 
lambdaL1H = p.Results.lambdaL1H; 
W_fixed = p.Results.W_fixed; 
W_init = p.Results.W_init; 
H_init = p.Results.H_init; 
SortFactors = p.Results.SortFactors; 
lambdaOrthoW = p.Results.lambdaOrthoW; 
lambdaOrthoH = p.Results.lambdaOrthoH; 
useWupdate = p.Results.useWupdate; 
M = p.Results.M; 

%test input
if isempty(spikes)
    error(['No spikes data']);
end

% working on this
% if ~isempty(Intervals)
%     t1=Intervals(1,1)*spikes.sr;
%     t2=Intervals(1,2)*spikes.sr;
%     bined_time = [t1:spikes.sr*bin:t2];
% else
    for i=1:length(spikes.ts)
        max_t(1,i) = max(spikes.ts{1,i});
    end  
    t1=0;t2=max(max_t);
    bined_time = [t1:spikes.sr*bin:t2];   
% end

%prepare spike train
for i=1:length(spikes.ts)
    spikes_tmp{1,i} = spikes.ts{1,i}(find(spikes.ts{1,i}>=t1 & spikes.ts{1,i}<=t2)); 
    X(i,:) = histc(spikes_tmp{1,i},bined_time);
end    

% Check that we have non-negative data
if min(X(:)) < 0
    error('Negative values in data!');
end

% zeropad data by L
X = [zeros(N,L),X,zeros(N,L)];
[N, T] = size(X);

%% initialize

% initialize W_init, H_init and M, if not provided
if isnan(W_init)
    W_init = max(X(:))*rand(N, K, L);
end
if isnan(H_init)
    H_init = max(X(:))*rand(K,T)./(sqrt(T/3)); % normalize so frobenius norm of each row ~ 1
% else
%     H_init = [zeros(K,L),H_init,zeros(K,L)];
end
if isnan(M)
    M = ones(N,T);
else
    M = [ones(N,L),M,ones(N,L)];
end

Xhat = helper.reconstruct(W, H); 
mask = find(M == 0); % find masked (held-out) indices 
X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit

smoothkernel = ones(1,(2*L)-1);  % for factor competition
smallnum = max(X(:))*1e-6; 
lasttime = 0;

% Calculate initial cost
cost = zeros(maxiter+1, 1);
cost(1) = sqrt(mean((X(:)-Xhat(:)).^2));

%% run for maxiter times
for iter = 1 : maxiter
    % Stopping criteria... Stop if reach maxiter or if change in cost function is less than the tolerance
    if (iter == maxiter) || ((iter>5) && (cost(iter+1)+tolerance)>mean(cost((iter-5):iter)))
        cost = cost(1 : iter+1);  % trim vector
        lasttime = 1; 
        if iter>1
            lambda = 0; % Do one final CNMF iteration (no regularization, just prioritize reconstruction)
        end
    end
    
    % Compute terms for standard CNMF H update 
    WTX = zeros(K, T);
    WTXhat = zeros(K, T);
    for l = 1 : L
        X_shifted = circshift(X,[0,-l+1]); 
        Xhat_shifted = circshift(Xhat,[0,-l+1]); 
        WTX = WTX + W(:, :, l)' * X_shifted;
        WTXhat = WTXhat + W(:, :, l)' * Xhat_shifted;
    end   
         
    % Compute regularization terms for H update
    if lambda>0
        dRdH = lambda.*(~eye(K))*conv2(WTX, smoothkernel, 'same');  
    else 
        dRdH = 0; 
    end
    if lambdaOrthoH>0
        dHHdH = lambdaOrthoH*(~eye(K))*conv2(H, smoothkernel, 'same');
    else
        dHHdH = 0;
    end
    dRdH = dRdH + lambdaL1H + dHHdH; % include L1 sparsity, if specified
    
    % Update H
    H = H .* WTX ./ (WTXhat + dRdH +eps);
        
    % Shift to center factors
    if shift
        [W, H] = helper.shiftFactors(W, H);  
        W = W+smallnum; % add small number to shifted W's, since multiplicative update cannot effect 0's
    end
    
    % Renormalize so rows of H have constant energy
    norms = sqrt(sum(H.^2, 2))';
    H = diag(1 ./ (norms+eps)) * H;
    for l = 1 : L
        W(:, :, l) = W(:, :, l) * diag(norms);
    end 
    
    if ~W_fixed
    % Update each Wl separately
        Xhat = helper.reconstruct(W, H); 
        mask = find(M == 0); % find masked (held-out) indices 
        X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit
        if lambdaOrthoW>0
            Wflat = sum(W,3);
        end
        if lambda>0 && useWupdate
            XS = conv2(X, smoothkernel, 'same'); 
        end
        for l = 1 : L % could parallelize to speed up for long L
            % Compute terms for standard CNMF W update
            H_shifted = circshift(H,[0,l-1]);
            XHT = X * H_shifted';
            XhatHT = Xhat * H_shifted';

            % Compute regularization terms for W update
            if lambda>0 && useWupdate; % Often get similar results with just H update, so option to skip W update
                dRdW = lambda.*XS*(H_shifted')*(~eye(K)); 
            else
                dRdW = 0;
            end
            if lambdaOrthoW>0
                dWWdW = lambdaOrthoW*Wflat*(~eye(K));
            else
                dWWdW = 0;
            end
            dRdW = dRdW + lambdaL1W + dWWdW; % include L1 and Worthogonality sparsity, if specified
            % Update W
            W(:, :, l) = W(:, :, l) .* XHT ./ (XhatHT + dRdW + eps);
        end
    end
    % Calculate cost for this iteration
    Xhat = helper.reconstruct(W, H);    
    mask = find(M == 0); % find masked (held-out) indices 
    X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit
    cost(iter+1) = sqrt(mean((X(:)-Xhat(:)).^2));

    % Plot to show progress
    if showPlot 
        SimpleWHPlot(W, H, Xhat,0); 
        title(sprintf('iteration #%i',iter));
        drawnow
    end
    
    if lasttime
        break
    end
end

% Undo zeropadding by truncating X, Xhat and H
X = X(:,L+1:end-L);
Xhat = Xhat(:,L+1:end-L);
H = H(:,L+1:end-L);

% Compute explained power of whole reconstruction and each factor
power = (sum(X(:).^2)-sum((X(:)-Xhat(:)).^2))/sum(X(:).^2);  % fraction power explained by whole reconstruction
[loadings,ind] = sort(helper.computeLoadingPercentPower(X,W,H),'descend'); % fraction power explained by each factor

% sort factors by loading power
if params.SortFactors
    W = W(:,ind,:);
    H = H(ind,:);
end
    
end
