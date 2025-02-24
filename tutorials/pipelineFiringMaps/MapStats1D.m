function stats = MapStats1D(map,varargin)

%MapStats1D - Compute statistics and find place fields from 1D linearized firing maps.
% 
%  USAGE
%
%    stats = MapStats1D(map)
%
%
%   INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'map'   map obtained using <a href="matlab:help Map">Map</a>
%    
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.15)
%     'minSize'     fields smaller than this percentage of the maze size 
%                   are considered spurious and ignored (default = 0.05)
%     'maxSize'     fields larger than this percentage of the maze size 
%                   are considered noise and ignored (default = 0.50)
%     'sepEdge'     fields with maximum Firing Rate closer to the edges less
%                   than this percentage of the maze size are ignored
%                   (default = 0.0)
%                   are considered noise and ignored (default = 0.50)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1 Hz)
%     'minPeak2nd'  for secondary place fields, peaks smaller than this 
%                   percentage of maximum Firing Rate along the maze are
%                   considered spurious and ignored (default 0.60)
%     'doPlot'   	Make and save plots (default: true) 
%    =========================================================================
%
%   OUTPUTS
%    stats.x             abscissa of the maximum value (in bins)
%    stats.y             ordinate of the maximum value (in bins)
%    stats.peak          in-field maximum value
%    stats.mean          in-field mean value
%    stats.size          field size (in bins)
%    stats.field         field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX        field x boundaries (in bins)
%    stats.fieldY        field y boundaries (in bins)
%    stats.specificity   spatial specificity (in bits, see Skaggs et al., 1993)
%
%    currently not dealing with circular data 

% Original code from MapStats.m by Michaël Zugaro, findPlaceFieldsAvg1D.m
% Antonio FR, 10/2019, and PlceCellInfo.m by Azahara Oliva
% Modified by Wenbo Tang, Feb 01, 20023

% Parse inputs 
p=inputParser;
addParameter(p,'threshold',0.15,@isnumeric);
addParameter(p,'minSize',0.05,@isnumeric);
addParameter(p,'maxSize',0.60,@isnumeric);
addParameter(p,'minPeak',2,@isnumeric);
addParameter(p,'minPeak2nd',0.6,@isnumeric);
addParameter(p,'sepEdge',0.0,@isnumeric);
addParameter(p,'doPlot', true, @islogical);

parse(p,varargin{:});
sizeMaze = length(map.z);
threshold = p.Results.threshold;
minSize = p.Results.minSize * sizeMaze;
maxSize = p.Results.maxSize * sizeMaze;
sepEdge = p.Results.sepEdge * sizeMaze;
minPeak = p.Results.minPeak;
minPeak2nd = p.Results.minPeak2nd;
doPlot = p.Results.doPlot;

%% Find place fields
% Default values
stats.x = NaN;
stats.field = [];
stats.size = 0;
stats.peak = 0;
stats.mean = 0;
stats.fieldX = [NaN NaN];
stats.specificity = nan;
stats.informationPerSpike = nan;
stats.informationPerSec = nan;
stats.sparsity = nan;
stats.selectivity = nan;

map_original = map;
x = 1:length(map.z);

% Compute the spatial specificity of the map, based on the formula proposed by Skaggs et al. (1993):
%  specificity = SUM { p(i) . lambda(i)/lambda . log2(lambda(i)/lambda) }
% Compute the spatial information, sparsity and selectivity
T = sum(map.time(:));
p_i = map.time/(T+eps); % Probability of the animal occupying bin 'i'
lambda_i = map.z;
lambda = lambda_i(:)'*p_i(:);

meanFiringRate = sum(sum(map.z.*map.time))./T;
logArg = map.z./meanFiringRate;

if T > 0 && lambda ~= 0
    stats.specificity = sum(sum(p_i.*lambda_i/lambda.*log2(lambda_i/lambda)));
    logArg(logArg == 0) = 1;
    stats.informationPerSpike  = sum(sum(p_i.*logArg.*log2(logArg))); % bits per spike.
    stats.informationPerSec = sum(sum(p_i.*map.z.*log2(logArg))); % bits per second.
    stats.sparsity = ((sum(sum(p_i.*map.z))).^2)/sum(sum(p_i.*(map.z.^2)));
    stats.selectivity = max(max(map.z))./meanFiringRate;
end


        
% Maximum FR along maze
maxFR = max(max(map.z));

% If there is no firing rate, go to next unit
if isempty(maxFR) || maxFR == 0 
    stats.field = logical(zeros(size(map.z)));
end

nBinsX = max([1 length(x)]);	% minimum number of bins is 1
circX = 0; circY = 0;

% Each time we find a field, we will remove it from the map; make a copy first
% Try to find more fields until no remaining bin exceeds min value
i = 1;
while true
    % Are there any candidate (unvisited) peaks left?
    [peak,idx] = max(map.z(:));
    % If separation from edges is less than sepEdge, go to next unit
    if (idx < sepEdge) || (idx > sizeMaze-sepEdge)
        break;
    end
    % If FR peak of 1st PF is less than minPeak, go to next unit
    % If FR peak of 2nd PF is less than minPeak2nd of maximum FR,
    % go to next unit
    if peak < ((i==1)*minPeak + (i>1)*maxFR*minPeak2nd)
    	break;
    end
    % Determine coordinates of largest candidate peak
    [y,x] = ind2sub(size(map.z),idx);
    % Find field (using min threshold for inclusion)
    field1 = FindFieldHelper(map.z,x,y,peak*threshold,circX,circY);
    size1 = sum(field1(:));
    % Does this field include two coalescent subfields?
    % To answer this question, we simply re-run the same field-searching procedure on the field
    % we then either keep the original field or choose the subfield if the latter is less than
    % 1/2 the size of the former
    m = peak*threshold;
    field2 = FindFieldHelper(map.z-m,x,y,(peak-m)*threshold,circX,circY);
    size2 = sum(field2(:));
    if size2 < 1/2*size1
        field = field2;
        tc = ' ';sc = '*'; % for debugging messages
    else
        field = field1;
        tc = '*';sc = ' '; % for debugging messages
    end
            
    % If rate map between place fields doesn't go below threshold,
    % discard new place field
    good2ndPF = true;
    if i > 1
        field0ini = find(diff(isnan(map.z))==1); if length(field0ini)>1, field0ini = field0ini(2); end
        field0end = find(diff(isnan(map.z))==-1); if length(field0end)>1, field0end = field0end(2); end
        field1ini = find(diff(field)==1); if isempty(field1ini), field1ini = 1; end
        field1end = find(diff(field)==-1);
        [~,idxBetwFields] = min([abs(field1ini-field0end),abs(field0ini-field1end)]);
        if idxBetwFields == 1
            if ~any(map.z(field1end:field0ini)<peak*threshold), good2ndPF = false; end
        else
            if ~any(map.z(field0end:field1ini)<peak*threshold), good2ndPF = false; end
        end
    end
            
    fieldSize = sum(field(:));
    % Keep this field if its size is sufficient
    if (fieldSize > minSize) && (fieldSize < maxSize) && good2ndPF
           stats.field(:,i) = field;
           stats.size(i) = fieldSize;
           stats.peak(i) = peak;
           stats.mean(i) = mean(map.z(field));
           idx = find(field & map.z == peak);
           [stats.y(i),stats.x(i)] = ind2sub(size(map.z),idx(1));
           [x,y] = FieldBoundaries(field,circX,circY);
           [stats.fieldX(i,:),stats.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
     end
     i = i + 1;
            
     % Mark field bins as visited
     map.z(field) = NaN;
     if all(isnan(map.z)), break; end
end

%%
% ==========
%   PLOT    
% ==========
if doPlot
    figure(1);
    plot(map_original.z,'k')
    if sum(map_original.z)>0
        hold on
        for ii = 1:size(stats.field,2)
            plot(find(stats.field(:,ii)),map_original.z(stats.field(:,ii) == 1),'linewidth',2)
            plot([1 1]*stats.x(ii),[0 map_original.z(stats.x(ii) == 1)],'--k')
        end
    end
    ylabel('FR(Hz)')
    xlabel('Track (#bin)')
    pause(0.4)
    clf
end
end
%%
% ------------------------------- Helper functions -------------------------------

% Field boundaries (circumscribed rectangle)

function [x,y] = FieldBoundaries(field,circX,circY)

% Find boundaries
x = find(any(field,1));
if isempty(x),
	x = [NaN NaN];
else
	x = [x(1) x(end)];
end
y = find(any(field,2));
if isempty(y),
	y = [NaN NaN];
else
	y = [y(1) y(end)];
end

% The above works in almost all cases; it fails however for circular coordinates if the field extends
% around an edge, e.g. for angles between 350° and 30°

if circX && x(1) == 1 && x(2) == size(field,2),
	xx = find(~all(field,1));
	if ~isempty(xx),
		x = [xx(end) xx(1)];
	end
end
if circY && y(1) == 1 && y(2) == size(field,1),
	yy = find(~all(field,2));
	if ~isempty(yy),
		y = [yy(end) yy(1)];
	end
end
end
