function [unitRip] = unitSWRmetrics(ripSpk,spikes,varargin)
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
%       FRevent         firing rate within each event
%       FRavgEvent      average firing rate outside events
%       
%
%
%   Antonio FR, 8/2020

p = inputParser;
addOptional(p,'baseFR',[])
parse(p,varargin{:})
baseFR = p.Results.baseFR;

%% Prob of participation
% for each unit
temp = zeros(length(ripSpk.UnitAbs),length(ripSpk.EventAbs));
for unit = 1:length(ripSpk.UnitAbs)
    for rip = 1:length(ripSpk.EventAbs)
        temp(unit,rip) = ~isempty(ripSpk.UnitEventAbs{unit,rip});
    end
end
unitRip.particip = sum(temp,2)/size(temp,2);

% for each event (n cells)
unitRip.nCellsEvent = (sum(temp,1))';
% get participating UIDs
nCellsEventUID = cell(size(temp,2),1);
temp = logical(temp);
for i = 1:size(temp,2)
    nCellsEventUID{i} = spikes.UID(temp(:,i));
end
unitRip.nCellsEventUID = nCellsEventUID;    
clear temp;
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
nUnitEvt = zeros(length(ripSpk.EventAbs),1);
unitUID = cell(length(ripSpk.EventAbs),1);
for unit = 1:length(ripSpk.UnitAbs)
    for rip = 1:length(ripSpk.EventAbs)
        unitRip.nSpkEach(unit,rip) = length(ripSpk.UnitEventAbs{unit,rip});
        unitTemp = unitUID{rip};
        if ~unitRip.nSpkEach(unit,rip)
            nUnitEvt(rip) = nUnitEvt(rip)+1;
            unitTemp(end+1) = spikes.UID(unit);
        end 
        unitUID{rip} = unitTemp;
    end
end
unitRip.nUnitEvt = nUnitEvt;    
unitRip.unitUID = unitUID;
   
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
for rip = 1:size(unitRip.nSpkEach,2)
    tempSum = [];
    for un = 1:size(unitRip.nSpkEach,1)
        if unitRip.nSpkEach(un,rip) ~= 0
            tempSum = [tempSum unitRip.nSpkEach(un,rip)];
        end
    end
    unitRip.FRevent(rip) = (sum(tempSum)/length(tempSum))/ripSpk.EventDuration(rip);
end