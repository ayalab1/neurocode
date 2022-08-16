function smoothed = smooth(data,smooth,varargin)

%Smooth - Smooth using a Gaussian kernel.
%
%  USAGE
%
%    smoothed = Smooth(data,smooth,<options>)
%
%    data           data to smooth
%    smooth         vertical and horizontal kernel size:
%                    - gaussian: standard deviations [Sv Sh]
%                    - rect/triangular: window half size [Nv Nh]
%                   (in number of samples, 0 = no smoothing)
%    <options>      optional list of property-value pairs (see table below))
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   - for 1D data, only one letter is used (default 'll')
%  	'kernel'		  either 'gaussian' (default), 'rect' (running average), 
%  					  or 'triangular' (weighted running average)
%    =========================================================================
%

% Copyright (C) 2004-2015 by Michaël Zugaro, 2013 Nicolas Maingret
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

maxSize = 10001;

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

vector = isvector(data);

if vector && any(isnan(data)),
    smoothed = data;
    nans = isnan(data);
    data(nans) = [];
    smoothed(~nans) = Smooth(data, smooth, varargin{:});
    return
end
matrix = (~vector & length(size(data)) == 2);
if ~vector & ~matrix,
	error('Smoothing applies only to vectors or matrices (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

% Vectors must be 'vertical'
if size(data,1) == 1,
	data = data';
end

% Default values
kernel = 'gaussian';
if vector, type = 'l'; else type = 'll'; end

if ~isdvector(smooth,'>=0') | (matrix & length(smooth) > 2) | (vector & length(smooth) ~= 1),
	error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

% If Sh = Sv = 0, no smoothing required
if all(smooth==0),
	smoothed = data;
	return
end

if length(smooth) == 1,
	% For 2D data, providing only one value S for the std is interpreted as Sh = Sv = S
	smooth = [smooth smooth];
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			type = lower(varargin{i+1});
			if (vector && ~isastring(type,'c','l')) || (~vector && ~isastring(type,'cc','cl','lc','ll')),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
			end
		case 'kernel',
			kernel = lower(varargin{i+1});
			if ~isastring(kernel,'gaussian','rectangular','triangular'),
				error('Incorrect value for property ''kernel'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).']);

  end
end

% Build kernels
[vSize,hSize] = size(data);

if strcmp(kernel,'rectangular'),

	% Rectangular kernel (running average)
	% 1) Vertical kernel
	vKernelSize = 2*smooth(1)+1;
	vKernel = ones(vKernelSize,1);
	vKernel = vKernel/sum(vKernel);
	% 2) Horizontal kernel
	if ~vector,
		hKernelSize = 2*smooth(2)+1;
		hKernel = ones(hKernelSize,1);
		hKernel = hKernel/sum(hKernel);
	end

elseif strcmp(kernel,'triangular'),

	% Triangular kernel
	% 1) Vertical kernel
	vKernelSize = 2*smooth(1)+1;
	vKernel = triang(vKernelSize);
	vKernel = vKernel/sum(vKernel);
	% 2) Horizontal kernel
	if ~vector,
		hKernelSize = 2*smooth(2)+1;
		hKernel = ones(hKernelSize,1);
		hKernel = hKernel/sum(hKernel);
	end

else

	% Gaussian kernel
	% 1) Vertical kernel
% 	if vSize > maxSize, warning(['Kernel too large; using ' int2str(maxSize) ' points.']); end
	vKernelSize = min([vSize maxSize]);
	r = (-vKernelSize:vKernelSize)'/vKernelSize;
	vKernelStdev = smooth(1)/vKernelSize;
	vKernel = exp(-r.^2/(vKernelStdev+eps)^2/2);
	vKernel = vKernel/sum(vKernel);
	% 2) Horizontal kernel
	if ~vector,
% 		if hSize > maxSize, warning(['Kernel too large; using ' int2str(maxSize) ' points.']); end
		hKernelSize = min([hSize maxSize]);
		r = (-hKernelSize:hKernelSize)/hKernelSize;
		hKernelStdev = smooth(2)/hKernelSize;
		hKernel = exp(-r.^2/(hKernelStdev+eps)^2/2);
		hKernel = hKernel/sum(hKernel);
	end
	
end

if vector,
	% Vector smoothing
	% Prepend/append data to limit edge effects
	if strcmp(type,'l'),
		% For linear data, flip edge data
		% top = 2*data(1)-flipud(data(1:vKernelSize));
		% bottom = 2*data(end)-flipud(data(end-vKernelSize+1:end));
		top = flipud(data(1:vKernelSize));
		bottom = flipud(data(end-vKernelSize+1:end));
	else
		% For circular data, wrap edge data
		top = data(end-vKernelSize+1:end);
		bottom = data(1:vKernelSize);
	end
	data = [top;data;bottom];
	% Convolve (and return central part)
	tmp = conv(vKernel,data);
	n = size(tmp,1);
	d = n - vSize;
	start = d/2+1;
	stop = start + vSize - 1;
	smoothed = tmp(start:stop,:);
else
	% Matrix smoothing
	% Convolve
	if smooth(1) == 0,    
		% Smooth only across columns (Sv = 0)
		% Prepend/append data to limit edge effects
		if strcmp(type(1),'l'),
            if any(isnan(data(:))),
                % take care of edges first
                indices = meshgrid(1:size(data,2),1:size(data,1));
                indices(isnan(data))  = nan;
                okLeft = min(indices,[],2);
                okRight = max(indices,[],2);
                for i=fliplr(2:max(okLeft)),
                    data(okLeft>=i,i-1) = data(okLeft>=i,i);
                end
                for i=(min(okRight):size(data,2)-1),
                    data(okRight<=i,i+1) = data(okRight<=i,i);
                end
                % interpolate mid-nans
                data = data';
                nans = isnan(data(:));
                indices = (1:numel(data))';
                data(nans) = interp1(indices(~nans),data(~nans),indices(nans));
                data = data';
            end
			% For linear data, flip edge data
			% left = 2*repmat(data(:,1),1,hKernelSize)-fliplr(data(:,1:hKernelSize));
			% right = 2*repmat(data(:,end),1,hKernelSize)-fliplr(data(:,end-hKernelSize+1:end));
			left = fliplr(data(:,1:hKernelSize));
			right = fliplr(data(:,end-hKernelSize+1:end));
        else
            if any(isnan(data(:))),
                % take care of edges first
                data = data';
                okLeft = find(~isnan(data(:)),1,'first');
                data(1:okLeft) = data(okLeft);
                okRight = find(~isnan(data(:)),1,'last');
                data(okRight:end) = data(okRight);
                % interpolate mid-nans
                nans = isnan(data(:));
                indices = (1:numel(data))';
                data(nans) = interp1(indices(~nans),data(~nans),indices(nans));
                data = data';
            end
			% For circular data, wrap edge data
			left = data(:,end-hKernelSize+1:end);
			right = data(:,1:hKernelSize);
		end
		data = [left data right];
		for i = 1:size(data,1),
			tmp = conv(hKernel,data(i,:));
			n = size(tmp,2);
			d = n - hSize;
			start = d/2+1;
			stop = start + hSize - 1;
			smoothed(i,:) = tmp(:,start:stop);
		end
	elseif smooth(2) == 0,
		% Smooth only across lines (Sh = 0)
		% Prepend/append data to limit edge effects
        if any(isnan(data(:))),
                % take care of edges first
                [~,indices] = meshgrid(1:size(data,2),1:size(data,1));
                indices(isnan(data))  = nan;
                okTop = min(indices);
                okBottom = max(indices);
                for i=fliplr(2:max(okTop)),
                    data(i-1,okTop>=i) = data(i,okTop>=i);
                end
                for i=(min(okBottom):size(data,1)-1),
                    data(i+1,okBottom<=i) = data(i,okBottom<=i);
                end
                % interpolate mid-nans
                nans = isnan(data(:));
                indices = (1:numel(data))';
                data(nans) = interp1(indices(~nans),data(~nans),indices(nans));
            end
		if strcmp(type(2),'l'),
			% For linear data, flip edge data
			% top = 2*repmat(2*data(1,:),vKernelSize,1)-flipud(data(1:vKernelSize,:));
			% bottom = 2*repmat(2*data(end,:),vKernelSize,1)-flipud(data(end-vKernelSize+1:end,:));
			top = flipud(data(1:vKernelSize,:));
			bottom = flipud(data(end-vKernelSize+1:end,:));
        else
            if any(isnan(data(:))),
                % take care of edges first
                okLeft = find(~isnan(data(:)),1,'first');
                data(1:okLeft) = data(okLeft);
                okRight = find(~isnan(data(:)),1,'last');
                data(okRight:end) = data(okRight);
                % interpolate mid-nans
                nans = isnan(data(:));
                indices = (1:numel(data))';
                data(nans) = interp1(indices(~nans),data(~nans),indices(nans));
            end
			% For circular data, wrap edge data
			bottom = data(1:vKernelSize,:);
			top = data(end-vKernelSize+1:end,:);
		end
		data = [top;data;bottom];
		for i = 1:size(data,2),
			tmp = conv(vKernel,data(:,i));
			n = size(tmp,1);
			d = n - vSize;
			start = d/2+1;
			stop = start + vSize - 1;
			smoothed(:,i) = tmp(start:stop);
		end
	else
        % Smooth in 2D
        if any(isnan(data(:))),
            if strcmp(type(1),'c') || strcmp(type(2),'c'),
                warning('linear interpolation for missing values');
            end
            % Note (raly): I know this step should be in the linear condition only,
            % but for the moment I haven' t had time to generalise the fillnans \
            % function to do circular interpolation
            data = fillnans(data,smooth);
        end
		% Prepend/append data to limit edge effects
		if strcmp(type(2),'l'),
			% For linear data, flip edge data
  			% top = 2*repmat(2*data(1,:),vKernelSize,1)-flipud(data(1:vKernelSize,:));
  			% bottom = 2*repmat(2*data(end,:),vKernelSize,1)-flipud(data(end-vKernelSize+1:end,:));
			top = flipud(data(1:vKernelSize,:));
			bottom = flipud(data(end-vKernelSize+1:end,:));
		else
			% For circular data, wrap edge data
			bottom = data(1:vKernelSize,:);
			top = data(end-vKernelSize+1:end,:);
		end
		data = [top;data;bottom];
		if strcmp(type(1),'l'),
			% For linear data, flip edge data
			% left = 2*repmat(2*data(:,1),1,hKernelSize)-fliplr(data(:,1:hKernelSize));
  			% right = 2*repmat(2*data(:,end),1,hKernelSize)-fliplr(data(:,end-hKernelSize+1:end));
			left = fliplr(data(:,1:hKernelSize));
			right = fliplr(data(:,end-hKernelSize+1:end));
		else
			% For circular data, wrap edge data
			left = data(:,end-hKernelSize+1:end);
			right = data(:,1:hKernelSize);
		end
		data = [left data right];
		tmp = conv2(vKernel,hKernel,data,'same');
		n = size(tmp,1);
		d = n - vSize;
		vStart = d/2+1;
		vStop = vStart + vSize - 1;
		n = size(tmp,2);
		d = n - hSize;
		hStart = d/2+1;
		hStop = hStart + hSize - 1;
		smoothed = tmp(vStart:vStop,hStart:hStop);
	end
end
end


%% ==  HELPER FUNCTION  ==

function interpolated = fillnans(matrix,coefficients)

% coeffifcients are vertical and horizontal kernel size (like in Smooth):
% how much should each dimension weigh in the interpolation?

% This function is an usage of John D'Errico's method 2 in
% INPAINT_NANS (dowloaded from FileExchange on 01/02/2016)
% e-mail address: woodchips@rochester.rr.com
% which can do exactly the same interpolation process for
% coefficients  =  [1 1];
% INPAINT_NANS: in-paints over nans in an array
% usage: B = INPAINT_NANS(A)          % default method
% usage: B = INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%         This method uses a simple plate metaphor.
%         del^2 is used.
%         Extrapolation behavior is linear.
%
%         Note: This method has problems in 1-d, so this
%         method is disabled for vector inputs.


% I always need to know which elements are NaN,
% and what size the array is for any method
[nRows,nCols] = size(matrix);
matrix = matrix(:);
nm = nRows*nCols;
k = isnan(matrix(:));

if nargin<2,
    coefficients = [1 1];
end

% list the nodes which are known, and which will
% be interpolated
nan_list = find(k);
known_list = find(~k);

% convert NaN indices to (r,c) form
% nan_list =  = find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc] = ind2sub([nRows,nCols],nan_list);

% both forms of index in one array:
% column 1  =  =  unrolled index
% column 2  =  =  row index
% column 3  =  =  column index
nan_list = [nan_list,nr,nc];

% Direct solve for del^2 BVP across holes

% generate sparse array with second partials on row
% variable for each nan element, only for those nodes
% which have a row index > 1 or < n


% a 2-d case
L  =  find((nan_list(:,2) > 1) & (nan_list(:,2) < nRows));
nl = length(L);
if nl>0
    fda = sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1]*coefficients(1),nl,1),nRows*nCols,nRows*nCols);
else
    fda = spalloc(nRows*nCols,n*nCols,size(nan_list,1)*5);
end

% 2nd partials on column index
L  =  find((nan_list(:,3) > 1) & (nan_list(:,3) < nCols));
nl = length(L);
if nl>0
    fda = fda+sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-nRows 0 nRows],nl,1), ...
        repmat([1 -2 1]*coefficients(2),nl,1),nRows*nCols,nRows*nCols);
end

% fix boundary conditions at extreme corners
% of the array in case there were nans there
if ismember(1,nan_list(:,1)),fda(1,[1 2 nRows+1]) = [-2 1 1];end
if ismember(nRows,nan_list(:,1)),fda(nRows,[nRows, nRows-1,nRows+nRows]) = [-2 1 1];end
if ismember(nm-nRows+1,nan_list(:,1)),fda(nm-nRows+1,[nm-nRows+1,nm-nRows+2,nm-nRows]) = [-2 1 1];end
if ismember(nm,nan_list(:,1)),fda(nm,[nm,nm-1,nm-nRows]) = [-2 1 1];end

% eliminate knowns
rhs = -fda(:,known_list)*matrix(known_list);

% and solve...
interpolated = matrix;
k = nan_list(:,1);
interpolated(k) = fda(k,k)\rhs(k);


% all done, make sure that B is the same shape as
% A was when we came in.
interpolated = reshape(interpolated,nRows,nCols);
end