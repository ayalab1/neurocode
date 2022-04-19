function simulated = GeneratePoissonSpikes(spikes,maxTime)

% Simulate poisson spikes with the same firing rate at the provided spikes
% spikes should be provided in [timestamp id] format

simulated = [];
for i=1:max(spikes(:,2))
    s = spikes(spikes(:,2)==i);
    isi = exprnd(maxTime/length(s),length(s),1);
    sh = cumsum(isi);
    sh(:,2) = i;
    simulated = [simulated;sh];
end
simulated(simulated(:,1)>maxTime,:) = [];
simulated = sortrows(simulated);
