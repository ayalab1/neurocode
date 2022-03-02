function [Information_1,Information_2,Sparsity,Coefficient,Selectivity] = PlaceCellInfo(Rate_Map, snSpikes, sTimeSpent)

% Compute different informative index about Place Cells
% INPUTS:
% 1- TimeSpent = Time Spent (s) per bin
% 2- nSpikes = Spikes (count) per bin
% 3- Rate_Map = nSpikes / TimeSpent
% OUTPUTS:
% 1- Information Index 1 = bits per spike
% 2- Information Index 2 = bits per second
% 3- Sparsity = 
% 4- Coefficient = index to quantify the sparsity (I am using directly the
%                  sparsity not this one)
% 5- Selectivity


if any(size(Rate_Map)~=size(snSpikes)) || any(size(Rate_Map)~=size(sTimeSpent))
    error('any(size(Rate_Map)~=size(snSpikes)) || any(size(Rate_Map)~=size(sTimeSpent))')
end

Information_1 = NaN;
Information_2 = NaN;
Sparsity = NaN;
Coefficient = NaN;
Selectivity = NaN;

T = sum(sum(sTimeSpent));
%meanFiringRate = sum(sum(map.firings))/(sum(sum(map.time)+eps));
if T > 0 & sum(sum(snSpikes))>0   
    occupancy = sTimeSpent./T;
    %meanFiringRate = sum(sum(snSpikes))./T;
    meanFiringRate = sum(sum(Rate_Map.*sTimeSpent))./T;
    
    logArg = Rate_Map./meanFiringRate;
%%% Is the following definition OK?
    logArg(logArg == 0) = 1;
        %%stats.information = sum(sum(snSpikes.*log2(logArg).*occupancy))/meanFiringRate;
    Information_1 = sum(sum(occupancy.*logArg.*log2(logArg))); % bits per spike.
    Information_2 = sum(sum(occupancy.*Rate_Map.*log2(logArg))); % bits per second.
        
    Sparsity = ((sum(sum(occupancy.*Rate_Map))).^2)/sum(sum(occupancy.*(Rate_Map.^2)));
        
    Coefficient = sqrt((1/Sparsity)-1);
    Selectivity = max(max(Rate_Map))./meanFiringRate;
end