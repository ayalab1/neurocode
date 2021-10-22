function [unitRip] = unitSWRmetrics(ripSpk,varargin)
%        [unitRip] = unitSWRmetrics(ripSpk)
%
%   INPUTS
%
%   ripSpk      buzcode ripple structure from getRipSpikes
%   baseFR      vector with baseline firing rates per unit. Optional, to caculate ripple gain
%
%   OUTPUTS
%
%
%
%   Antonio FR, 8/2020

p = inputParser;
addParameter(p,'baseFR',[])
parse(p,varargin{:})
baseFR = p.Results.baseFR;

%% Prob of participation
% for each unit
for unit = 1:length(ripSpk.UnitAbs)
    for rip = 1:length(ripSpk.EventAbs)
        if ~isempty(ripSpk.UnitEventAbs{unit,rip})
            temp(unit,rip) = 1;
        else
            temp(unit,rip) = 0;
        end
    end
end
unitRip.particip = sum(temp,2)/size(temp,2);

% for each event (n cells)
unitRip.nCellsEvent = (sum(temp,1)/size(temp,1))'; clear temp;

%% Firing rates
% in each ripple
for unit = 1:length(ripSpk.UnitAbs)
    for rip = 1:length(ripSpk.EventAbs)
        unitRip.FReach(unit,rip) = length(ripSpk.UnitEventAbs{unit,rip})/ripSpk.EventDuration(rip);
    end
end
% mean all ripples
for unit = 1:length(ripSpk.UnitAbs)
    unitRip.FRall(unit,1) = length(ripSpk.UnitAbs{unit})/sum(ripSpk.EventDuration);
end
% population mean all ripples
for rip = 1:length(ripSpk.EventAbs)
    tmp=0;
    for unit = 1:length(ripSpk.UnitAbs)
        tmp=tmp+length(ripSpk.UnitEventAbs{unit,rip});
    end
    unitRip.FRpopulation(rip,1) = tmp/sum(ripSpk.EventDuration);
end
% mean all ripples in which unit fired
for unit = 1:length(ripSpk.UnitAbs)
    unitRip.FRparticip(unit,1) = mean(nonzeros(unitRip.FReach(unit,:)));
end

%% Gain
if baseFR
    % in each ripple
    for unit = 1:length(ripSpk.UnitAbs)
        for rip = 1:length(ripSpk.EventAbs)
            unitRip.GainEach(unit,rip) = unitRip.FReach(unit,rip)/baseFR(unit,1);
        end
    end
    % mean all ripples
    for unit = 1:length(ripSpk.UnitAbs)
        unitRip.GainAll(unit,1) = mean(unitRip.GainEach(unit,:));
    end
    % mean all ripples in which unit fired
    for unit = 1:length(ripSpk.UnitAbs)
        unitRip.GainParticip(unit,1) = mean(nonzeros(unitRip.GainEach(unit,:)));
    end
end

%% n spikes
% in each ripple
for unit = 1:length(ripSpk.UnitAbs)
    for rip = 1:length(ripSpk.EventAbs)
        unitRip.nSpkEach(unit,rip) = length(ripSpk.UnitEventAbs{unit,rip});
    end
end
% mean all ripples
for unit = 1:length(ripSpk.UnitAbs)
    unitRip.nSpkAll(unit,1) = mean((unitRip.nSpkEach(unit,:)));
end
% mean all ripples in which unit fired
for unit = 1:length(ripSpk.UnitAbs)
    unitRip.nSpkParticip(unit,1) = mean(nonzeros(unitRip.nSpkEach(unit,:)));
end
% mean all units each event
unitRip.nSpkEvent = sum(unitRip.nSpkEach,1)';

% mean FR per event
for rip = 1:length(ripSpk.EventAbs)
    unitRip.FRevent = unitRip.nSpkEvent/ripSpk.EventDuration(rip);
end

end

