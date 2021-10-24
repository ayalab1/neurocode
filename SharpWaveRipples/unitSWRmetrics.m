function [unitRip] = unitSWRmetrics(ripSpk,varargin)
%        [unitRip] = unitSWRmetrics(ripSpk)
%
%   INPUTS
%
%   ripSpk      buzcode ripple structure from getRipSpikes
%   baseFR      vector with baseline firing rates per unit. Optional, to caculate ripple gain
%
%   OUTPUTS
%   unitRip
%       particip        fraction of events that each unit participates in
%       nCellsEvent     fraction of cells that participate in each event
%       FReach          firing rate for each unit in each event
%       FRall           firing rate for each unit during events
%       FRparticip      firing rate during events where the cell fired
%       nSpkEach        number of spikes in each unit and each event
%       nSpkAll         number of spikes for each unit during events
%       nSpkParticip    number of spikes during events where the cell fired
%       nSpkEvent       number of spikes during each event
%       nSpkNonEvent    number of spikes outside each event
%       FRevent         firing rate within each event
%       FRnonEvent      firing rate outside each event
%       FRavgNon        average firing rate outside events
%       FRavgEvent      average firing rate outside events
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
% mean all ripples in which unit fired
for unit = 1:length(ripSpk.UnitAbs)
    unitRip.FRparticip(unit,1) = mean(nonzeros(unitRip.FReach(unit,:)));
end
% outside of events
for unit = 1:length(ripSpk.UnitNonAbs)
    for rip = 1:length(ripSpk.NonEventAbs)
        unitRip.FReachNon(unit,rip) = length(ripSpk.UnitNonEventAbs{unit,rip})/ripSpk.NonEventDur(rip);
    end
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

%% n spikes non
% outside each ripple
for unit = 1:length(ripSpk.UnitNonAbs)
    for rip = 1:length(ripSpk.NonEventAbs)
        unitRip.nSpkEachNon(unit,rip) = length(ripSpk.UnitNonEventAbs{unit,rip});
    end
end
% mean all ripples
for unit = 1:length(ripSpk.UnitNonAbs)
    unitRip.nSpkAllNon(unit,1) = mean((unitRip.nSpkEachNon(unit,:)));
end
% mean all ripples in which unit fired
for unit = 1:length(ripSpk.UnitNonAbs)
    unitRip.nSpkParticipNon(unit,1) = mean(nonzeros(unitRip.nSpkEachNon(unit,:)));
end
% mean all units each event
unitRip.nSpkNonEvent = sum(unitRip.nSpkEachNon,1)';

% mean FR per event
for rip = 1:length(ripSpk.NonEventAbs)
    unitRip.FRnonEvent = unitRip.nSpkNonEvent/ripSpk.NonEventDur(rip);
end
end