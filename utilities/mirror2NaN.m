function [outputvar] = mirror2NaN(inputvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
% inputvar = var of any type (with values of any kind)
%
%
% OUTPUTS
% outputvar = same class but = NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aza 2022

%check emptyness
if isempty(inputvar)
    outputvar=[];
    warning('Just so you know your input was empty and so will be your output');
    return
end

%main stuff, select the type of data
if isstruct(inputvar)
    
    outputvar = nanstruct(inputvar);
    return
elseif iscell(inputvar)
    
    outputvar = nancell(inputvar);
    return
else
    
    outputvar = nanmat(inputvar);
    return
end

%individual functions for nan different types of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%struct class
function outstruct = nanstruct(instruct)
    names = fieldnames(instruct);
    
    for i=1:length(names)
        if isstruct(eval(['instruct.' names{i} ]))
            
           outstruct.(eval('names{i}')) = nanstruct(['instruct.' names{i} ]);
        
        elseif iscell(eval(['instruct.' names{i} ]))
            
           outstruct.(eval('names{i}')) = nancell(eval(['instruct.' names{i} ])); 
                        
        else
            
           outstruct.(eval('names{i}')) = nanmat(eval(['instruct.' names{i} ])); 
            
        end
    end
end

%cell class
function outcell = nancell(incell)
    for dim1 =1:size(incell,1)
        for dim2 =1:size(incell,2)
            outcell{dim1,dim2} = NaN*incell{dim1,dim2};
        end
    end
end

%mat (any, in theory...) class
function outmat = nanmat(inmat)
	outmat = NaN*inmat;
end

end