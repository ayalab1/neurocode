function [ structout ] = collapseStruct( structin,dim,combine,NEST )
%structout = CollapseStruct( structin,dim,combine,NEST ) Combines elements in a
%structure array
%
%INPUT
%   structin    struct(N).fields structure array with N elements where each
%               of the N elements has the same fields and field structure.
%               structin can have nested cell arrays or stuctures, but be
%               careful with this... may not work for >1-dimensional
%               structures/cell arrays
%   dim         dimension along which to combine each element in the
%               structure array. (default: 2)
%               option: 'match'. CollapseStruct will find the dimension
%               along which the elements are concatenatable. If there are
%               multiple, will pick the first one.
%   (optional)
%       combine     can take 'mean' or 'median' instead of concatenating
%                   or 'std'. default: 'justcat'
%       NEST        true/false - go into nested structures?
%
%       
%% A place for input options
if ~exist('NEST','var')
    NEST = false;
end

if ~exist('combine','var')
    combine = 'justcat';
end

if ~exist('dim','var')
    dim = 2;
end

%%
if ~isstruct(structin)
    structout=structin;
    warning('Not a structure...')
    return
end
fields = fieldnames(structin);

%%

for ff = 1:length(fields)
    currentfield = fields{ff};

	if isstruct(structin(1).(currentfield)) & NEST %For Nested Structures
        %Here: need to match fields, either 'remove' or 'addempty'
        %bz_Matchfields
       %structout.(currentfield) = cat(1,structin(:).(currentfield));
       structout.(currentfield) = matchfields({structin(:).(currentfield)},[],'remove');
       try
       structout.(currentfield) = collapseStruct(structout.(currentfield),dim,combine,true);
       catch
           keyboard
       end
       continue
    elseif iscell(structin(1).(currentfield)) & NEST %For cell array in field
        if strcmp(dim,'match')
            catdim = bz_FindCatableDims({structin(:).(currentfield)});
            if length(catdim)>1
                catdim = catdim(1); %Change this to be the largest dimension
            end
        else
            catdim = dim;
        end
        try
            structout.(currentfield) = cat(catdim,structin(:).(currentfield));
        catch
            display(['Failed to concatenate field ',currentfield])
            %keyboard
        end
    elseif (isstring(structin(1).(currentfield))||ischar(structin(1).(currentfield))) & NEST %For string in field
        structout.(currentfield) = {structin(:).(currentfield)};
    else %For simple array in field
        try
            %Concatenate Array
            if strcmp(dim,'match')
                catdim = findCatableDims({structin(:).(currentfield)});
            else
                catdim = dim;
            end
        structout.(currentfield) = cat(catdim(1),structin(:).(currentfield));
        catch
             display(['Failed to concatenate field ',currentfield])
           % keyboard
%             continue
         end
	end

    
        if ~exist('structout','var')
            structout = [];
            continue
        elseif ~isfield(structout,currentfield)
            continue
        end
    switch combine
        case 'mean'
            structout.(currentfield) = nanmean(structout.(currentfield),catdim);
        case 'median'
            structout.(currentfield) = nanmedian(structout.(currentfield),catdim);
        case 'std'
            structout.(currentfield) = nanstd(structout.(currentfield),[],catdim);
        case 'sem'
            structout.(currentfield) =  nanstd(structout.(currentfield),[],catdim)./sqrt(sum(~isnan(structout.(currentfield)),catdim));
    end
end

function [ catabledim ] = findCatableDims( tryarrays )
%Remove any empty arrays
tryarrays(cellfun(@(X) isempty(X),tryarrays)) = [];

arraysizes = cellfun(@size,tryarrays,'uniformoutput',false);
numdims = cellfun(@length,arraysizes);

%If all the arrays don't have the same number of dimensions, you're out of luck
if range(numdims)~=0
    disp('Your arrays don''t have the same number of dimensions. You''re SOL')
end

arraysizes = cat(1,arraysizes{:});

samesizes = range(arraysizes,1)==0;

if all(samesizes)
    catabledim = 1:numdims;
elseif sum(~samesizes)==1
    catabledim = find(~samesizes);
else
    disp('Two or more of your dimensions are different size')
end

end

end

