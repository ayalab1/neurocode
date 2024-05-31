function [templates,correlations,weights,variance] = ActivityTemplates(spikes,varargin)

%ActivityTemplates - Compute assemblies activity and neurons contribution from PCA of spike trains and ICA.
%
% Computes the templates for the component activation analysis described in
% Peyrache et al (2009). Optionally, computes ICA to improve detection. These
% templates can then be tested on different data sets using <a href="matlab:help ReactivationStrength">ReactivationStrength</a>.
% Time bins can be automatically determined using a fixed bin size or provided
% as an explicit list (e.g. computed using theta phases).
%
%  USAGE
%
%    [templates,correlations,weights,nSignNeurons] = ActivityTemplates(spikes,<options>)
%
%    spikes         spike train (either single unit or MUA)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bins'        list of [start stop] for all bins
%     'binSize'     bin size in s (default = 0.050)
%     'step'		step size in s (default = same as binSize)
%     'mode'		algorithm to perform ('pca' or 'ica'; default = 'ica')
%     'tracyWidom'  whether the Tracy Widom correction should be used when
%                   selecting for significant PCs (default = false)
%    =========================================================================
%
%  OUTPUT
%
%    templates      3D array of template matrices (dimension 3 is template #,
%                   ordered in descending order of corresponding eigenvalue
%                   for 'pca')
%    correlations   a neuron-by-neuron correlation matrix over the provided
%                   bins, serving as basis for the analysis
%    weights		cell weights from for each assemblly; they correspond to 
%                   eigenvectors in 'pca' mode and ICs in 'ica' mode
%    variance       proportion of the variance explained by the weights; they
%                   correspond to eigenvalues in 'pca' mode.
%    
%  SEE
%
%    See also ReactivationStrength.

% Copyright (C) 2016-2022 by Michaël Zugaro, Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
bins = [];
defaultBinSize = 0.050;
binSize = [];
step = [];
mode = 'ica';
tracyWidom = false;
controlBins = [];

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
end
% Check parameter sizes
if ~isdmatrix(spikes,'@2'),
	error('Parameter ''spikes'' is not a Nx2 matrix (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
			end
		case 'step',
			step = varargin{i+1};
			if ~isdscalar(step,'>0'),
				error('Incorrect value for property ''step'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
			end
		case 'bins',
			bins = varargin{i+1};
			if ~isdmatrix(bins,'@2'),
				error('Incorrect value for property ''bins'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
			end
		case 'mode',
			mode = varargin{i+1};
			if ~isastring(mode,'pca','ica','varimax'),
				error('Incorrect value for property ''bins'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
            end
        case 'tracywidom',
            tracywidom = lower(varargin{i+1});
            % allow for 'on'/'off' usage and transform it to true/false
            if isastring(tracyWidom,'on','off'), if strcmp(tracyWidom,'on'), tracyWidom = true; else tracyWidom = false; end; end
            if length(tracyWidom)==1
                if ~islogical(tracyWidom), tracyWidom = tracyWidom>0; end % allow for 0/1 usage and transform it to true/false
            else
                error('Incorrect value for property ''bins'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
            end
        case 'controlbins',
            controlBins = varargin{i+1};
            if ~isdmatrix(controlBins,'@2'),
                error('Incorrect value for property ''controlBins'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).']);
    end
end

% Options binSize and bins are incompatible
if ~isempty(binSize) && ~isempty(bins),
	error('Parameters ''binSize'' and ''bins'' are incompatible (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
end
if isempty(binSize) && isempty(bins),
	binSize = defaultBinSize;
end
if isempty(step), step = binSize; end
nUnits = max(spikes(:,2));
templates = nan(nUnits,nUnits,0);
correlations = nan(nUnits,nUnits);
weights = nan(nUnits,0);
if isempty(nUnits), return; end

%% Bin spikes
spikes = sortrows(spikes,1);
id = spikes(:,2);

% Shift spike times to start at 0, and list bins unless explicitly provided
if isempty(bins),
	spikes(:,1) = spikes(:,1) - spikes(1,1);
	bins = (0:step:(spikes(end,1)-binSize))';
	bins(:,2) = bins+binSize;
else
	m = min([min(spikes(:,1)) min(bins(:))]);
	spikes(:,1) = spikes(:,1) - m;
	bins = bins - m;
end

% Create spike count matrix
nBins = size(bins,1);
if isempty(nBins), return; end
n = zeros(nBins,nUnits);
for unit = 1:nUnits,
	n(:,unit) = CountInIntervals(spikes(id==unit,1),bins);
end

% Crate a control spike count matrix constructing correlations to ignore:
if ~isempty(controlBins)
    n0 = zeros(size(controlBins,1),nUnits);
    for unit = 1:nUnits,
    	n0(:,unit) = CountInIntervals(spikes(id==unit,1),controlBins);
    end
    controlCorrelations = (1/(nBins-1))*n0'*n0;
end

%% Create correlation matrix
n = zscore(n);
correlations = (1/(nBins-1))*n'*n;

if ~isempty(controlBins)
    correlations = correlations - controlCorrelations;
end

% Compute eigenvalues/vectors and sort in descending order
[eigenvectors,eigenvalues] = eig(correlations);
[eigenvalues,i] = sort(diag(eigenvalues),'descend');
eigenvectors = eigenvectors(:,i);
eigenvectors(:,max(eigenvectors)~=max(abs(eigenvectors))) = -eigenvectors(:,max(eigenvectors)~=max(abs(eigenvectors))); % flip so that largest weight is positive

%% Keep only significant eigenvalues and compute templates

q = nBins/nUnits;
if q < 1,
	warning('Not enough time bins to determine significant templates');
	eigenvalues = NaN;
end

lambdaMax = (1+sqrt(1/q))^2;
if tracyWidom
    lambdaMax =  lambdaMax + nUnits^(-2/3); % Tracy-Widom correction
end
significant = eigenvalues>lambdaMax;

if sum(significant)==0,
	templates = zeros(nUnits,nUnits,0); weights = zeros(nUnits,0); return;
end

eigenvectors = eigenvectors(:,significant);

if strcmp(mode,'pca'),
    templates = zeros(nUnits,nUnits,sum(significant));
	weights = eigenvectors;
	for i = 1:sum(significant),
		templates(:,:,i) = weights(:,i)*weights(:,i)';
		templates(:,:,i) = templates(:,:,i) - diag(diag(templates(:,:,i))); % remove the diagonal
    end
    variance = eigenvalues(significant)/nUnits;
	return
end

if strcmp(mode,'varimax'),
    try
    [weights,~] = rotatefactors(eigenvectors(:,significant),'method','varimax');
    catch
        warning('Varimax did not work. Returning eigevectors straight out of PCA.');
        weights = eigenvectors(:,significant);
    end
    templates = zeros(nUnits,nUnits,sum(significant));
%     The sign of the weights in a component is arbitrary (+component and -component are equivalent)
%     Flip weights so that the most deviating weight of a component is positive (it is more convenient for visualisation that the assembly has positive weights)
    flip = max(weights)<-min(weights);
    weights(:,flip) = -weights(:,flip);
    for i = 1:size(weights,2),
        templates(:,:,i) = weights(:,i)*weights(:,i)';
        templates(:,:,i) = templates(:,:,i) - diag(diag(templates(:,:,i))); % remove the diagonal
    end
    variance = var(n*weights)/nUnits;
    return
end

projection = (eigenvectors * eigenvectors') * n';

%% Run the ICA on the new spike matrix
[weights,~] = fastica(projection,'pcaE',eigenvectors,'pcaD',diag(eigenvalues(significant)));

if isempty(weights),
	templates = zeros(nUnits,nUnits,0); weights = zeros(nUnits,0); return;
end

% The sign of the weights in a component is arbitrary (+component and -component are equivalent)
% Flip weights so that the most deviating weight of a component is positive (it is more convenient for visualisation that the assembly has positive weights)
flip = max(weights)<-min(weights);
weights(:,flip) = -weights(:,flip);

% Normalise weights as Van de Ven et al (2016):
nor = 1; for i=1:size(weights,2),nor(1,i)=norm(weights(:,i)); end
weights = bsxfun(@rdivide,weights,nor);


variance = var(n*weights)/nUnits; % the total variance of the activity matrix "n" is nUnits because the matrix is z-scored (variance of 1 per column)
% Order them by the amount of variance they explain (for consistency's sake):
[~,order] = sort(-variance); % from highest to lowest
weights = weights(:,order);
templates = zeros(nUnits,nUnits,size(weights,2));
for i = 1:size(weights,2)
	templates(:,:,i) = weights(:,i)*weights(:,i)';
	templates(:,:,i) = templates(:,:,i) - diag(diag(templates(:,:,i))); % remove the diagonal
end

% ------------------------------- Helper functions -------------------------------

function [Out1, Out2, Out3] = fastica(mixedsig, varargin)
%FASTICA - Fast Independent Component Analysis
%
% FastICA for Matlab 7.x and 6.x
% Version 2.5, October 19 2005
% Copyright (c) Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.
%
% FASTICA(mixedsig) estimates the independent components from given
% multidimensional signals. Each row of matrix mixedsig is one
% observed signal.  FASTICA uses Hyvarinen's fixed-point algorithm,
% see http://www.cis.hut.fi/projects/ica/fastica/. Output from the
% function depends on the number output arguments:
%
% [icasig] = FASTICA (mixedsig); the rows of icasig contain the
% estimated independent components.
%
% [icasig, A, W] = FASTICA (mixedsig); outputs the estimated separating
% matrix W and the corresponding mixing matrix A.
%
% [A, W] = FASTICA (mixedsig); gives only the estimated mixing matrix
% A and the separating matrix W.
%
% Some optional arguments induce other output formats, see below.
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
% FASTICA can be called with numerous optional arguments. Optional
% arguments are given in parameter pairs, so that first argument is
% the name of the parameter and the next argument is the value for
% that parameter. Optional parameter pairs can be given in any order.
%
% OPTIONAL PARAMETERS:
%
% Parameter name        Values and description
%
%======================================================================
% --Basic parameters in fixed-point algorithm:
%
% 'approach'            (string) The decorrelation approach used. Can be
%                       symmetric ('symm'), i.e. estimate all the
%                       independent component in parallel, or
%                       deflation ('defl'), i.e. estimate independent
%                       component one-by-one like in projection pursuit.
%                       Default is 'defl'.
%
% 'numOfIC'             (integer) Number of independent components to
%                       be estimated. Default equals the dimension of data.
%
%======================================================================
% --Choosing the nonlinearity:
%
% 'g'                   (string) Chooses the nonlinearity g used in
%                       the fixed-point algorithm. Possible values:
%
%                       Value of 'g':      Nonlinearity used:
%                       'pow3' (default)   g(u)=u^3
%                       'tanh'             g(u)=tanh(a1*u)
%                       'gauss             g(u)=u*exp(-a2*u^2/2)
%                       'skew'             g(u)=u^2
%
% 'finetune'		(string) Chooses the nonlinearity g used when
%                       fine-tuning. In addition to same values
%                       as for 'g', the possible value 'finetune' is:
%                       'off'              fine-tuning is disabled.
%
% 'a1'                  (number) Parameter a1 used when g='tanh'.
%                       Default is 1.
% 'a2'                  (number) Parameter a2 used when g='gaus'.
%                       Default is 1.
%
% 'mu'			(number) Step size. Default is 1.
%                       If the value of mu is other than 1, then the
%                       program will use the stabilized version of the
%                       algorithm (see also parameter 'stabilization').
%
%
% 'stabilization'       (string) Values 'on' or 'off'. Default 'off'.
%                       This parameter controls wether the program uses
%                       the stabilized version of the algorithm or
%                       not. If the stabilization is on, then the value
%                       of mu can momentarily be halved if the program
%                       senses that the algorithm is stuck between two
%                       points (this is called a stroke). Also if there
%                       is no convergence before half of the maximum
%                       number of iterations has been reached then mu
%                       will be halved for the rest of the rounds.
%
%======================================================================
% --Controlling convergence:
%
% 'epsilon'             (number) Stopping criterion. Default is 0.0001.
%
% 'maxNumIterations'    (integer) Maximum number of iterations.
%                       Default is 1000.
%
% 'maxFinetune'         (integer) Maximum number of iterations in
%                       fine-tuning. Default 100.
%
% 'sampleSize'          (number) [0 - 1] Percentage of samples used in
%                       one iteration. Samples are chosen in random.
%                       Default is 1 (all samples).
%
% 'initGuess'           (matrix) Initial guess for A. Default is random.
%                       You can now do a "one more" like this:
%                       [ica, A, W] = fastica(mix, 'numOfIC',3);
%                       [ica2, A2, W2] = fastica(mix, 'initGuess', A, 'numOfIC', 4);
%
%======================================================================
% --Graphics and text output:
%
% 'verbose'             (string) Either 'on' or 'off'. Default is
%                       'on': report progress of algorithm in text format.
%
% 'displayMode'         (string) Plot running estimates of independent
%                       components: 'signals', 'basis', 'filters' or
%                       'off'. Default is 'off'.
%
% 'displayInterval'     Number of iterations between plots.
%                       Default is 1 (plot after every iteration).
%
%======================================================================
% --Controlling reduction of dimension and whitening:
%
% Reduction of dimension is controlled by 'firstEig' and 'lastEig', or
% alternatively by 'interactivePCA'.
%
% 'firstEig'            (integer) This and 'lastEig' specify the range for
%                       eigenvalues that are retained, 'firstEig' is
%                       the index of largest eigenvalue to be
%                       retained. Default is 1.
%
% 'lastEig'             (integer) This is the index of the last (smallest)
%                       eigenvalue to be retained. Default equals the
%                       dimension of data.
%
% 'interactivePCA'      (string) Either 'on' or 'off'. When set 'on', the
%                       eigenvalues are shown to the user and the
%                       range can be specified interactively. Default
%                       is 'off'. Can also be set to 'gui'. Then the user
%                       can use the same GUI that's in FASTICAG.
%
% If you already know the eigenvalue decomposition of the covariance
% matrix, you can avoid computing it again by giving it with the
% following options:
%
% 'pcaE'                (matrix) Eigenvectors
% 'pcaD'                (matrix) Eigenvalues
%
% If you already know the whitened data, you can give it directly to
% the algorithm using the following options:
%
% 'whiteSig'            (matrix) Whitened signal
% 'whiteMat'            (matrix) Whitening matrix
% 'dewhiteMat'          (matrix) dewhitening matrix
%
% If values for all the 'whiteSig', 'whiteSig' and 'dewhiteMat' are
% supplied, they will be used in computing the ICA. PCA and whitening
% are not performed. Though 'mixedsig' is not used in the main
% algorithm it still must be entered - some values are still
% calculated from it.
%
% Performing preprocessing only is possible by the option:
%
% 'only'                (string) Compute only PCA i.e. reduction of
%                       dimension ('pca') or only PCA plus whitening
%                       ('white'). Default is 'all': do ICA estimation
%                       as well.  This option changes the output
%                       format accordingly. For example:
%
%                       [whitesig, WM, DWM] = FASTICA(mixedsig,
%                       'only', 'white')
%                       returns the whitened signals, the whitening matrix
%                       (WM) and the dewhitening matrix (DWM). (See also
%                       WHITENV.) In FastICA the whitening matrix performs
%                       whitening and the reduction of dimension. Dewhitening
%                       matrix is the pseudoinverse of whitening matrix.
%
%                       [E, D] = FASTICA(mixedsig, 'only', 'pca')
%                       returns the eigenvector (E) and diagonal
%                       eigenvalue (D) matrices  containing the
%                       selected subspaces.
%
%======================================================================
% EXAMPLES
%
%       [icasig] = FASTICA (mixedsig, 'approach', 'symm', 'g', 'tanh');
%               Do ICA with tanh nonlinearity and in parallel (like
%               maximum likelihood estimation for supergaussian data).
%
%       [icasig] = FASTICA (mixedsig, 'lastEig', 10, 'numOfIC', 3);
%               Reduce dimension to 10, and estimate only 3
%               independent components.
%
%       [icasig] = FASTICA (mixedsig, 'verbose', 'off', 'displayMode', 'off');
%               Don't output convergence reports and don't plot
%               independent components.
%
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
%   See also FASTICAG

% @(#)$Id: fastica.m,v 1.14 2005/10/19 13:05:34 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check some basic requirements of the data
if nargin == 0,
	error ('You must supply the mixed data as input argument.');
end

if length (size (mixedsig)) > 2,
	error ('Input data can not have more than two dimensions.');
end

if any (any (isnan (mixedsig))),
	error ('Input data contains NaN''s.');
end

if ~isa (mixedsig, 'double')
	fprintf ('Warning: converting input data into regular (double) precision.\n');
	mixedsig = double (mixedsig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data

mixedmean = mean(mixedsig,2);
mixedsig = bsxfun(@minus,mixedsig,mixedmean);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values for optional parameters

% All
verbose           = 'off';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'off';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for fastICA - i.e. this file

b_verbose = 0;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters

if (rem(length(varargin),2)==1)
	error('Optional parameters should always go by pairs');
else
	for i=1:2:(length(varargin)-1)
		if ~ischar (varargin{i}),
			error (['Unknown type of optional parameter name (parameter' ...
				' names must be strings).']);
		end
		% change the value of parameter
		switch lower (varargin{i})
			case 'stabilization'
				stabilization = lower (varargin{i+1});
			case 'maxfinetune'
				maxFinetune = varargin{i+1};
			case 'samplesize'
				sampleSize = varargin{i+1};
			case 'verbose'
				verbose = lower (varargin{i+1});
				% silence this program also
				if strcmp (verbose, 'off'), b_verbose = 0; end
			case 'firsteig'
				firstEig = varargin{i+1};
			case 'lasteig'
				lastEig = varargin{i+1};
			case 'interactivepca'
				interactivePCA = lower (varargin{i+1});
			case 'approach'
				approach = lower (varargin{i+1});
			case 'numofic'
				numOfIC = varargin{i+1};
				% User has supplied new value for numOfIC.
				% We'll use this information later on...
				userNumOfIC = 1;
			case 'g'
				g = lower (varargin{i+1});
			case 'finetune'
				finetune = lower (varargin{i+1});
			case 'a1'
				a1 = varargin{i+1};
			case 'a2'
				a2 = varargin{i+1};
			case {'mu', 'myy'}
				myy = varargin{i+1};
			case 'epsilon'
				epsilon = varargin{i+1};
			case 'maxnumiterations'
				maxNumIterations = varargin{i+1};
			case 'initguess'
				% no use setting 'guess' if the 'initState' is not set
				initState = 'guess';
				guess = varargin{i+1};
			case 'displaymode'
				displayMode = lower (varargin{i+1});
			case 'displayinterval'
				displayInterval = varargin{i+1};
			case 'pcae'
				% calculate if there are enought parameters to skip PCA
				jumpPCA = jumpPCA + 1;
				E = varargin{i+1};
			case 'pcad'
				% calculate if there are enought parameters to skip PCA
				jumpPCA = jumpPCA + 1;
				D = varargin{i+1};
			case 'whitesig'
				% calculate if there are enought parameters to skip PCA and whitening
				jumpWhitening = jumpWhitening + 1;
				whitesig = varargin{i+1};
			case 'whitemat'
				% calculate if there are enought parameters to skip PCA and whitening
				jumpWhitening = jumpWhitening + 1;
				whiteningMatrix = varargin{i+1};
			case 'dewhitemat'
				% calculate if there are enought parameters to skip PCA and whitening
				jumpWhitening = jumpWhitening + 1;
				dewhiteningMatrix = varargin{i+1};
			case 'only'
				% if the user only wants to calculate PCA or...
				switch lower (varargin{i+1})
					case 'pca'
						only = 1;
					case 'white'
						only = 2;
					case 'all'
						only = 3;
				end
				
			otherwise
				% Hmmm, something wrong with the parameter string
				error(['Unrecognized parameter: ''' varargin{i} '''']);
		end;
	end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print information about data
if b_verbose
	fprintf('Number of signals: %d\n', Dim);
	fprintf('Number of samples: %d\n', NumOfSampl);
end

% Check if the data has been entered the wrong way,
% but warn only... it may be on purpose

if Dim > NumOfSampl
	if b_verbose
		fprintf('Warning: ');
		fprintf('The signal matrix may be oriented in the wrong way.\n');
		fprintf('In that case transpose the matrix.\n\n');
		input('Proceed anyway? (any key but Ctrl-C)');
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PCA

% We need the results of PCA for whitening, but if we don't
% need to do whitening... then we dont need PCA...
if jumpWhitening == 3
	if b_verbose,
		fprintf ('Whitened signal and corresponding matrises supplied.\n');
		fprintf ('PCA calculations not needed.\n');
	end;
else
	
	% OK, so first we need to calculate PCA
	% Check to see if we already have the PCA data
	if jumpPCA == 2,
		if b_verbose,
			fprintf ('Values for PCA calculations supplied.\n');
			fprintf ('PCA calculations not needed.\n');
		end;
	else
		% display notice if the user entered one, but not both, of E and D.
		if (jumpPCA > 0) & (b_verbose),
			fprintf ('You must suply all of these in order to jump PCA:\n');
			fprintf ('''pcaE'', ''pcaD''.\n');
		end;
		
		% Calculate PCA
		[E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
	end
end

% skip the rest if user only wanted PCA
if only > 1
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Whitening the data
	
	% Check to see if the whitening is needed...
	if jumpWhitening == 3,
		if b_verbose,
			fprintf ('Whitening not needed.\n');
		end;
	else
		
		% Whitening is needed
		% display notice if the user entered some of the whitening info, but not all.
		if (jumpWhitening > 0) & (b_verbose),
			fprintf ('You must suply all of these in order to jump whitening:\n');
			fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
		end;
		
		% Calculate the whitening
		whiteningMatrix = inv(sqrt(D))*E';
		dewhiteningMatrix = E*sqrt(D);
		whitesig =  whiteningMatrix*mixedsig;
		
	end
	
end % if only > 1

% skip the rest if user only wanted PCA and whitening
if only > 2
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculating the ICA
	
	% Check some parameters
	% The dimension of the data may have been reduced during PCA calculations.
	% The original dimension is calculated from the data by default, and the
	% number of IC is by default set to equal that dimension.
	
	Dim = size(whitesig, 1);
	
	% The number of IC's must be less or equal to the dimension of data
	if numOfIC > Dim
		numOfIC = Dim;
		% Show warning only if verbose = 'on' and user supplied a value for 'numOfIC'
		if (b_verbose & userNumOfIC)
			fprintf('Warning: estimating only %d independent components\n', numOfIC);
			fprintf('(Can''t estimate more independent components than dimension of data)\n');
		end
	end
	
	%% Calculate the ICA with fixed point algorithm.
	[A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
		numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
		maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
		displayMode, displayInterval, verbose);
	
	% Check for valid return
	if ~isempty(W)
		% Add the mean back in.
		if b_verbose
			fprintf('Adding the mean back to the data.\n');
		end
		icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl);
		%icasig = W * mixedsig;
		if b_verbose & ...
				(max(abs(W * mixedmean)) > 1e-9) & ...
				(strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
			fprintf('Note that the plots don''t have the mean added.\n');
		end
	else
		icasig = [];
	end
	
end % if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output depends on the number of output parameters
% and the 'only' parameter.

if only == 1    % only PCA
	Out1 = E;
	Out2 = D;
elseif only == 2  % only PCA & whitening
	if nargout == 2
		Out1 = whiteningMatrix;
		Out2 = dewhiteningMatrix;
	else
		Out1 = whitesig;
		Out2 = whiteningMatrix;
		Out3 = dewhiteningMatrix;
	end
else      % ICA
	if nargout == 2
		Out1 = A;
		Out2 = W;
	else
		Out1 = icasig;
		Out2 = A;
		Out3 = W;
	end
end


function [A, W] = fpica(X, whiteningMatrix, dewhiteningMatrix, approach, ...
	numOfIC, g, finetune, a1, a2, myy, stabilization, ...
	epsilon, maxNumIterations, maxFinetune, initState, ...
	guess, sampleSize, displayMode, displayInterval, ...
	s_verbose);
%FPICA - Fixed point ICA. Main algorithm of FASTICA.
%
% [A, W] = fpica(whitesig, whiteningMatrix, dewhiteningMatrix, approach,
%        numOfIC, g, finetune, a1, a2, mu, stabilization, epsilon,
%        maxNumIterations, maxFinetune, initState, guess, sampleSize,
%        displayMode, displayInterval, verbose);
%
% Perform independent component analysis using Hyvarinen's fixed point
% algorithm. Outputs an estimate of the mixing matrix A and its inverse W.
%
% whitesig                              :the whitened data as row vectors
% whiteningMatrix                       :can be obtained with function whitenv
% dewhiteningMatrix                     :can be obtained with function whitenv
% approach      [ 'symm' | 'defl' ]     :the approach used (deflation or symmetric)
% numOfIC       [ 0 - Dim of whitesig ] :number of independent components estimated
% g             [ 'pow3' | 'tanh' |     :the nonlinearity used
%                 'gaus' | 'skew' ]
% finetune      [same as g + 'off']     :the nonlinearity used in finetuning.
% a1                                    :parameter for tuning 'tanh'
% a2                                    :parameter for tuning 'gaus'
% mu                                    :step size in stabilized algorithm
% stabilization [ 'on' | 'off' ]        :if mu < 1 then automatically on
% epsilon                               :stopping criterion
% maxNumIterations                      :maximum number of iterations
% maxFinetune                           :maximum number of iteretions for finetuning
% initState     [ 'rand' | 'guess' ]    :initial guess or random initial state. See below
% guess                                 :initial guess for A. Ignored if initState = 'rand'
% sampleSize    [ 0 - 1 ]               :percentage of the samples used in one iteration
% displayMode   [ 'signals' | 'basis' | :plot running estimate
%                 'filters' | 'off' ]
% displayInterval                       :number of iterations we take between plots
% verbose       [ 'on' | 'off' ]        :report progress in text format
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%       [A, W] = fpica(nv, wm, dwm);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also FASTICA, FASTICAG, WHITENV, PCAMAT

% @(#)$Id: fpica.m,v 1.7 2005/06/16 12:52:55 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variable for stopping the ICA calculations from the GUI
global g_FastICA_interrupt;
if isempty(g_FastICA_interrupt)
	clear global g_FastICA_interrupt;
	interruptible = 0;
else
	interruptible = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values

if nargin < 3, error('Not enough arguments!'); end
[vectorSize, numSamples] = size(X);
if nargin < 20, s_verbose = 'on'; end
if nargin < 19, displayInterval = 1; end
if nargin < 18, displayMode = 'on'; end
if nargin < 17, sampleSize = 1; end
if nargin < 16, guess = 1; end
if nargin < 15, initState = 'rand'; end
if nargin < 14, maxFinetune = 100; end
if nargin < 13, maxNumIterations = 1000; end
if nargin < 12, epsilon = 0.0001; end
if nargin < 11, stabilization = 'on'; end
if nargin < 10, myy = 1; end
if nargin < 9, a2 = 1; end
if nargin < 8, a1 = 1; end
if nargin < 7, finetune = 'off'; end
if nargin < 6, g = 'pow3'; end
if nargin < 5, numOfIC = vectorSize; end     % vectorSize = Dim
if nargin < 4, approach = 'defl'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the data

if ~isreal(X)
	error('Input has an imaginary part.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for verbose

switch lower(s_verbose)
	case 'on'
		b_verbose = 1;
	case 'off'
		b_verbose = 0;
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for approach

switch lower(approach)
	case 'symm'
		approachMode = 1;
	case 'defl'
		approachMode = 2;
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''approach''\n', approach));
end
if b_verbose, fprintf('Used approach [ %s ].\n', approach); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for numOfIC

if vectorSize < numOfIC
	error('Must have numOfIC <= Dimension!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the sampleSize
if sampleSize > 1
	sampleSize = 1;
	if b_verbose
		fprintf('Warning: Setting ''sampleSize'' to 1.\n');
	end
elseif sampleSize < 1
	if (sampleSize * numSamples) < 1000
		sampleSize = min(1000/numSamples, 1);
		if b_verbose
			fprintf('Warning: Setting ''sampleSize'' to %0.3f (%d samples).\n', ...
				sampleSize, floor(sampleSize * numSamples));
		end
	end
end
if b_verbose
	if  b_verbose & (sampleSize < 1)
		fprintf('Using about %0.0f%% of the samples in random order in every step.\n',sampleSize*100);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for nonlinearity.

switch lower(g)
	case 'pow3'
		gOrig = 10;
	case 'tanh'
		gOrig = 20;
	case {'gaus', 'gauss'}
		gOrig = 30;
	case 'skew'
		gOrig = 40;
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''g''\n', g));
end
if sampleSize ~= 1
	gOrig = gOrig + 2;
end
if myy ~= 1
	gOrig = gOrig + 1;
end

if b_verbose,
	fprintf('Used nonlinearity [ %s ].\n', g);
end

finetuningEnabled = 1;
switch lower(finetune)
	case 'pow3'
		gFine = 10 + 1;
	case 'tanh'
		gFine = 20 + 1;
	case {'gaus', 'gauss'}
		gFine = 30 + 1;
	case 'skew'
		gFine = 40 + 1;
	case 'off'
		if myy ~= 1
			gFine = gOrig;
		else
			gFine = gOrig + 1;
		end
		finetuningEnabled = 0;
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''finetune''\n', ...
			finetune));
end

if b_verbose & finetuningEnabled
	fprintf('Finetuning enabled (nonlinearity: [ %s ]).\n', finetune);
end

switch lower(stabilization)
	case 'on'
		stabilizationEnabled = 1;
	case 'off'
		if myy ~= 1
			stabilizationEnabled = 1;
		else
			stabilizationEnabled = 0;
		end
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''stabilization''\n', ...
			stabilization));
end

if b_verbose & stabilizationEnabled
	fprintf('Using stabilized algorithm.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some other parameters
myyOrig = myy;
% When we start fine-tuning we'll set myy = myyK * myy
myyK = 0.01;
% How many times do we try for convergence until we give up.
failureLimit = 5;


usedNlinearity = gOrig;
stroke = 0;
notFine = 1;
long = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for initial state.

switch lower(initState)
	case 'rand'
		initialStateMode = 0;
	case 'guess'
		if size(guess,1) ~= size(whiteningMatrix,2)
			initialStateMode = 0;
			if b_verbose
				fprintf('Warning: size of initial guess is incorrect. Using random initial guess.\n');
			end
		else
			initialStateMode = 1;
			if size(guess,2) < numOfIC
				if b_verbose
					fprintf('Warning: initial guess only for first %d components. Using random initial guess for others.\n', size(guess,2));
				end
				guess(:, size(guess, 2) + 1:numOfIC) = ...
					rand(vectorSize,numOfIC-size(guess,2))-.5;
			elseif size(guess,2)>numOfIC
				guess=guess(:,1:numOfIC);
				fprintf('Warning: Initial guess too large. The excess column are dropped.\n');
			end
			if b_verbose, fprintf('Using initial guess.\n'); end
		end
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''initState''\n', initState));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for display mode.

switch lower(displayMode)
	case {'off', 'none'}
		usedDisplay = 0;
	case {'on', 'signals'}
		usedDisplay = 1;
		if (b_verbose & (numSamples > 10000))
			fprintf('Warning: Data vectors are very long. Plotting may take long time.\n');
		end
		if (b_verbose & (numOfIC > 25))
			fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
		end
	case 'basis'
		usedDisplay = 2;
		if (b_verbose & (numOfIC > 25))
			fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
		end
	case 'filters'
		usedDisplay = 3;
		if (b_verbose & (vectorSize > 25))
			fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
		end
	otherwise
		error(sprintf('Illegal value [ %s ] for parameter: ''displayMode''\n', displayMode));
end

% The displayInterval can't be less than 1...
if displayInterval < 1
	displayInterval = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_verbose, fprintf('Starting ICA calculation...\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYMMETRIC APPROACH
if approachMode == 1,
	
	% set some parameters more...
	usedNlinearity = gOrig;
	stroke = 0;
	notFine = 1;
	long = 0;
	
	A = zeros(vectorSize, numOfIC);  % Dewhitened basis vectors.
	if initialStateMode == 0
		% Take random orthonormal initial vectors.
		B = orth (randn (vectorSize, numOfIC));
	elseif initialStateMode == 1
		% Use the given initial vector as the initial state
		B = whiteningMatrix * guess;
	end
	
	BOld = zeros(size(B));
	BOld2 = zeros(size(B));
	
	% This is the actual fixed-point iteration loop.
	for round = 1:maxNumIterations + 1,
		if round == maxNumIterations + 1,
			fprintf('No convergence after %d steps\n', maxNumIterations);
			fprintf('Note that the plots are probably wrong.\n');
			if ~isempty(B)
				% Symmetric orthogonalization.
				B = B * real(inv(B' * B)^(1/2));
				
				W = B' * whiteningMatrix;
				A = dewhiteningMatrix * B;
			else
				W = [];
				A = [];
			end
			return;
		end
		
		if (interruptible & g_FastICA_interrupt)
			if b_verbose
				fprintf('\n\nCalculation interrupted by the user\n');
			end
			if ~isempty(B)
				W = B' * whiteningMatrix;
				A = dewhiteningMatrix * B;
			else
				W = [];
				A = [];
			end
			return;
		end
		
		
		% Symmetric orthogonalization.
		B = B * real(inv(B' * B)^(1/2));
		
		% Test for termination condition. Note that we consider opposite
		% directions here as well.
		minAbsCos = min(abs(diag(B' * BOld)));
		minAbsCos2 = min(abs(diag(B' * BOld2)));
		
		if (1 - minAbsCos < epsilon)
			if finetuningEnabled & notFine
				if b_verbose, fprintf('Initial convergence, fine-tuning: \n'); end;
				notFine = 0;
				usedNlinearity = gFine;
				myy = myyK * myyOrig;
				BOld = zeros(size(B));
				BOld2 = zeros(size(B));
				
			else
				if b_verbose, fprintf('Convergence after %d steps\n', round); end
				
				% Calculate the de-whitened vectors.
				A = dewhiteningMatrix * B;
				break;
			end
		elseif stabilizationEnabled
			if (~stroke) & (1 - minAbsCos2 < epsilon)
				if b_verbose, fprintf('Stroke!\n'); end;
				stroke = myy;
				myy = .5*myy;
				if mod(usedNlinearity,2) == 0
					usedNlinearity = usedNlinearity + 1;
				end
			elseif stroke
				myy = stroke;
				stroke = 0;
				if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
					usedNlinearity = usedNlinearity - 1;
				end
			elseif (~long) & (round>maxNumIterations/2)
				if b_verbose, fprintf('Taking long (reducing step size)\n'); end;
				long = 1;
				myy = .5*myy;
				if mod(usedNlinearity,2) == 0
					usedNlinearity = usedNlinearity + 1;
				end
			end
		end
		
		BOld2 = BOld;
		BOld = B;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Show the progress...
		if b_verbose
			if round == 1
				fprintf('Step no. %d\n', round);
			else
				fprintf('Step no. %d, change in value of estimate: %.3g \n', round, 1 - minAbsCos);
			end
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Also plot the current state...
		switch usedDisplay
			case 1
				if rem(round, displayInterval) == 0,
					% There was and may still be other displaymodes...
					% 1D signals
					icaplot('dispsig',(X'*B)');
					drawnow;
				end
			case 2
				if rem(round, displayInterval) == 0,
					% ... and now there are :-)
					% 1D basis
					A = dewhiteningMatrix * B;
					icaplot('dispsig',A');
					drawnow;
				end
			case 3
				if rem(round, displayInterval) == 0,
					% ... and now there are :-)
					% 1D filters
					W = B' * whiteningMatrix;
					icaplot('dispsig',W);
					drawnow;
				end
			otherwise
		end
		
		switch usedNlinearity
			% pow3
			case 10
				B = (X * (( X' * B) .^ 3)) / numSamples - 3 * B;
			case 11
				% optimoitu - epsilonin kokoisia eroja
				% t�m� on optimoitu koodi, katso vanha koodi esim.
				% aikaisemmista versioista kuten 2.0 beta3
				Y = X' * B;
				Gpow3 = Y .^ 3;
				Beta = sum(Y .* Gpow3);
				D = diag(1 ./ (Beta - 3 * numSamples));
				B = B + myy * B * (Y' * Gpow3 - diag(Beta)) * D;
			case 12
				Xsub=X(:, getSamples(numSamples, sampleSize));
				B = (Xsub * (( Xsub' * B) .^ 3)) / size(Xsub,2) - 3 * B;
			case 13
				% Optimoitu
				Ysub=X(:, getSamples(numSamples, sampleSize))' * B;
				Gpow3 = Ysub .^ 3;
				Beta = sum(Ysub .* Gpow3);
				D = diag(1 ./ (Beta - 3 * size(Ysub', 2)));
				B = B + myy * B * (Ysub' * Gpow3 - diag(Beta)) * D;
				
				% tanh
			case 20
				hypTan = tanh(a1 * X' * B);
				B = X * hypTan / numSamples - ...
					ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / numSamples * ...
					a1;
			case 21
				% optimoitu - epsilonin kokoisia
				Y = X' * B;
				hypTan = tanh(a1 * Y);
				Beta = sum(Y .* hypTan);
				D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
				B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
			case 22
				Xsub=X(:, getSamples(numSamples, sampleSize));
				hypTan = tanh(a1 * Xsub' * B);
				B = Xsub * hypTan / size(Xsub, 2) - ...
					ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / size(Xsub, 2) * a1;
			case 23
				% Optimoitu
				Y = X(:, getSamples(numSamples, sampleSize))' * B;
				hypTan = tanh(a1 * Y);
				Beta = sum(Y .* hypTan);
				D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
				B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
				
				% gauss
			case 30
				U = X' * B;
				Usquared=U .^ 2;
				ex = exp(-a2 * Usquared / 2);
				gauss =  U .* ex;
				dGauss = (1 - a2 * Usquared) .*ex;
				B = X * gauss / numSamples - ...
					ones(size(B,1),1) * sum(dGauss)...
					.* B / numSamples ;
			case 31
				% optimoitu
				Y = X' * B;
				ex = exp(-a2 * (Y .^ 2) / 2);
				gauss = Y .* ex;
				Beta = sum(Y .* gauss);
				D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex)));
				B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
			case 32
				Xsub=X(:, getSamples(numSamples, sampleSize));
				U = Xsub' * B;
				Usquared=U .^ 2;
				ex = exp(-a2 * Usquared / 2);
				gauss =  U .* ex;
				dGauss = (1 - a2 * Usquared) .*ex;
				B = Xsub * gauss / size(Xsub,2) - ...
					ones(size(B,1),1) * sum(dGauss)...
					.* B / size(Xsub,2) ;
			case 33
				% Optimoitu
				Y = X(:, getSamples(numSamples, sampleSize))' * B;
				ex = exp(-a2 * (Y .^ 2) / 2);
				gauss = Y .* ex;
				Beta = sum(Y .* gauss);
				D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex)));
				B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
				
				% skew
			case 40
				B = (X * ((X' * B) .^ 2)) / numSamples;
			case 41
				% Optimoitu
				Y = X' * B;
				Gskew = Y .^ 2;
				Beta = sum(Y .* Gskew);
				D = diag(1 ./ (Beta));
				B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;
			case 42
				Xsub=X(:, getSamples(numSamples, sampleSize));
				B = (Xsub * ((Xsub' * B) .^ 2)) / size(Xsub,2);
			case 43
				% Uusi optimoitu
				Y = X(:, getSamples(numSamples, sampleSize))' * B;
				Gskew = Y .^ 2;
				Beta = sum(Y .* Gskew);
				D = diag(1 ./ (Beta));
				B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;
				
			otherwise
				error('Code for desired nonlinearity not found!');
		end
	end
	
	
	% Calculate ICA filters.
	W = B' * whiteningMatrix;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Also plot the last one...
	switch usedDisplay
		case 1
			% There was and may still be other displaymodes...
			% 1D signals
			icaplot('dispsig',(X'*B)');
			drawnow;
		case 2
			% ... and now there are :-)
			% 1D basis
			icaplot('dispsig',A');
			drawnow;
		case 3
			% ... and now there are :-)
			% 1D filters
			icaplot('dispsig',W);
			drawnow;
		otherwise
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFLATION APPROACH
if approachMode == 2
	
	B = zeros(vectorSize);
	
	% The search for a basis vector is repeated numOfIC times.
	round = 1;
	
	numFailures = 0;
	
	while round <= numOfIC,
		myy = myyOrig;
		usedNlinearity = gOrig;
		stroke = 0;
		notFine = 1;
		long = 0;
		endFinetuning = 0;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Show the progress...
		if b_verbose, fprintf('IC %d ', round); end
		
		% Take a random initial vector of lenght 1 and orthogonalize it
		% with respect to the other vectors.
		if initialStateMode == 0
			w = randn (vectorSize, 1);
		elseif initialStateMode == 1
			w=whiteningMatrix*guess(:,round);
		end
		w = w - B * B' * w;
		w = w / norm(w);
		
		wOld = zeros(size(w));
		wOld2 = zeros(size(w));
		
		% This is the actual fixed-point iteration loop.
		%    for i = 1 : maxNumIterations + 1
		i = 1;
		gabba = 1;
		while i <= maxNumIterations + gabba
			if (usedDisplay > 0)
				drawnow;
			end
			if (interruptible & g_FastICA_interrupt)
				if b_verbose
					fprintf('\n\nCalculation interrupted by the user\n');
				end
				return;
			end
			
			% Project the vector into the space orthogonal to the space
			% spanned by the earlier found basis vectors. Note that we can do
			% the projection with matrix B, since the zero entries do not
			% contribute to the projection.
			w = w - B * B' * w;
			w = w / norm(w);
			
			if notFine
				if i == maxNumIterations + 1
					if b_verbose
						fprintf('\nComponent number %d did not converge in %d iterations.\n', round, maxNumIterations);
					end
					round = round - 1;
					numFailures = numFailures + 1;
					if numFailures > failureLimit
						if b_verbose
							fprintf('Too many failures to converge (%d). Giving up.\n', numFailures);
						end
						if round == 0
							A=[];
							W=[];
						end
						return;
					end
					% numFailures > failurelimit
					break;
				end
				% i == maxNumIterations + 1
			else
				% if notFine
				if i >= endFinetuning
					wOld = w; % So the algorithm will stop on the next test...
				end
			end
			% if notFine
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Show the progress...
			if b_verbose, fprintf('.'); end;
			
			
			% Test for termination condition. Note that the algorithm has
			% converged if the direction of w and wOld is the same, this
			% is why we test the two cases.
			if norm(w - wOld) < epsilon | norm(w + wOld) < epsilon
				if finetuningEnabled & notFine
					if b_verbose, fprintf('Initial convergence, fine-tuning: '); end;
					notFine = 0;
					gabba = maxFinetune;
					wOld = zeros(size(w));
					wOld2 = zeros(size(w));
					usedNlinearity = gFine;
					myy = myyK * myyOrig;
					
					endFinetuning = maxFinetune + i;
					
				else
					numFailures = 0;
					% Save the vector
					B(:, round) = w;
					
					% Calculate the de-whitened vector.
					A(:,round) = dewhiteningMatrix * w;
					% Calculate ICA filter.
					W(round,:) = w' * whiteningMatrix;
					
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% Show the progress...
					if b_verbose, fprintf('computed ( %d steps ) \n', i); end
					
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% Also plot the current state...
					switch usedDisplay
						case 1
							if rem(round, displayInterval) == 0,
								% There was and may still be other displaymodes...
								% 1D signals
								temp = X'*B;
								icaplot('dispsig',temp(:,1:numOfIC)');
								drawnow;
							end
						case 2
							if rem(round, displayInterval) == 0,
								% ... and now there are :-)
								% 1D basis
								icaplot('dispsig',A');
								drawnow;
							end
						case 3
							if rem(round, displayInterval) == 0,
								% ... and now there are :-)
								% 1D filters
								icaplot('dispsig',W);
								drawnow;
							end
					end
					% switch usedDisplay
					break; % IC ready - next...
				end
				%if finetuningEnabled & notFine
			elseif stabilizationEnabled
				if (~stroke) & (norm(w - wOld2) < epsilon | norm(w + wOld2) < ...
						epsilon)
					stroke = myy;
					if b_verbose, fprintf('Stroke!'); end;
					myy = .5*myy;
					if mod(usedNlinearity,2) == 0
						usedNlinearity = usedNlinearity + 1;
					end
				elseif stroke
					myy = stroke;
					stroke = 0;
					if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
						usedNlinearity = usedNlinearity - 1;
					end
				elseif (notFine) & (~long) & (i > maxNumIterations / 2)
					if b_verbose, fprintf('Taking long (reducing step size) '); end;
					long = 1;
					myy = .5*myy;
					if mod(usedNlinearity,2) == 0
						usedNlinearity = usedNlinearity + 1;
					end
				end
			end
			
			wOld2 = wOld;
			wOld = w;
			
			switch usedNlinearity
				% pow3
				case 10
					w = (X * ((X' * w) .^ 3)) / numSamples - 3 * w;
				case 11
					EXGpow3 = (X * ((X' * w) .^ 3)) / numSamples;
					Beta = w' * EXGpow3;
					w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
				case 12
					Xsub=X(:,getSamples(numSamples, sampleSize));
					w = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2) - 3 * w;
				case 13
					Xsub=X(:,getSamples(numSamples, sampleSize));
					EXGpow3 = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2);
					Beta = w' * EXGpow3;
					w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
					% tanh
				case 20
					hypTan = tanh(a1 * X' * w);
					w = (X * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / numSamples;
				case 21
					hypTan = tanh(a1 * X' * w);
					Beta = w' * X * hypTan;
					w = w - myy * ((X * hypTan - Beta * w) / ...
						(a1 * sum((1-hypTan .^2)') - Beta));
				case 22
					Xsub=X(:,getSamples(numSamples, sampleSize));
					hypTan = tanh(a1 * Xsub' * w);
					w = (Xsub * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / size(Xsub, 2);
				case 23
					Xsub=X(:,getSamples(numSamples, sampleSize));
					hypTan = tanh(a1 * Xsub' * w);
					Beta = w' * Xsub * hypTan;
					w = w - myy * ((Xsub * hypTan - Beta * w) / ...
						(a1 * sum((1-hypTan .^2)') - Beta));
					% gauss
				case 30
					% This has been split for performance reasons.
					u = X' * w;
					u2=u.^2;
					ex=exp(-a2 * u2/2);
					gauss =  u.*ex;
					dGauss = (1 - a2 * u2) .*ex;
					w = (X * gauss - sum(dGauss)' * w) / numSamples;
				case 31
					u = X' * w;
					u2=u.^2;
					ex=exp(-a2 * u2/2);
					gauss =  u.*ex;
					dGauss = (1 - a2 * u2) .*ex;
					Beta = w' * X * gauss;
					w = w - myy * ((X * gauss - Beta * w) / ...
						(sum(dGauss)' - Beta));
				case 32
					Xsub=X(:,getSamples(numSamples, sampleSize));
					u = Xsub' * w;
					u2=u.^2;
					ex=exp(-a2 * u2/2);
					gauss =  u.*ex;
					dGauss = (1 - a2 * u2) .*ex;
					w = (Xsub * gauss - sum(dGauss)' * w) / size(Xsub, 2);
				case 33
					Xsub=X(:,getSamples(numSamples, sampleSize));
					u = Xsub' * w;
					u2=u.^2;
					ex=exp(-a2 * u2/2);
					gauss =  u.*ex;
					dGauss = (1 - a2 * u2) .*ex;
					Beta = w' * Xsub * gauss;
					w = w - myy * ((Xsub * gauss - Beta * w) / ...
						(sum(dGauss)' - Beta));
					% skew
				case 40
					w = (X * ((X' * w) .^ 2)) / numSamples;
				case 41
					EXGskew = (X * ((X' * w) .^ 2)) / numSamples;
					Beta = w' * EXGskew;
					w = w - myy * (EXGskew - Beta*w)/(-Beta);
				case 42
					Xsub=X(:,getSamples(numSamples, sampleSize));
					w = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
				case 43
					Xsub=X(:,getSamples(numSamples, sampleSize));
					EXGskew = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
					Beta = w' * EXGskew;
					w = w - myy * (EXGskew - Beta*w)/(-Beta);
					
				otherwise
					error('Code for desired nonlinearity not found!');
			end
			
			% Normalize the new w.
			w = w / norm(w);
			i = i + 1;
		end
		round = round + 1;
	end
	if b_verbose, fprintf('Done.\n'); end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Also plot the ones that may not have been plotted.
	if (usedDisplay > 0) & (rem(round-1, displayInterval) ~= 0)
		switch usedDisplay
			case 1
				% There was and may still be other displaymodes...
				% 1D signals
				temp = X'*B;
				icaplot('dispsig',temp(:,1:numOfIC)');
				drawnow;
			case 2
				% ... and now there are :-)
				% 1D basis
				icaplot('dispsig',A');
				drawnow;
			case 3
				% ... and now there are :-)
				% 1D filters
				icaplot('dispsig',W);
				drawnow;
			otherwise
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the end let's check the data for some security
if ~isreal(A)
	if b_verbose, fprintf('Warning: removing the imaginary part from the result.\n'); end
	A = real(A);
	W = real(W);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction
% Calculates tanh simplier and faster than Matlab tanh.
function y=tanh(x)
y = 1 - 2 ./ (exp(2 * x) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Samples = getSamples(max, percentage)
Samples = find(rand(1, max) < percentage);


