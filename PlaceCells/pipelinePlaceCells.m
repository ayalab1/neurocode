%% Basic pipeline for place field analysis
% These are the basic fucntions for performing place field analysis (in
% 1D for now). Is still work in progress and anyone is welcome to contribute.

% TODO:
% - Decide a common trial structure and way to linearize positions.
% - Implent trial based firing map and place field detection
% - 2D place fields
% - phase precession
% - theta compression and sequences
% - position decoding

% 0- Linearize maze positions and get trial structure

% 1-Calculate firing maps
[firingMaps] = firingMapAvg(posTrials,spikes);
save([basename '.firingMapsAvg.cellinfo.mat'],'firingMaps');

% 1.1- save struct for cell explorer (you can display ratemaps in cell explorer)
ratemap = firingMaps;
for i = 1:length(ratemap.rateMaps)
    ratemap.map{i} = ratemap.rateMaps{i}{1};
end
ratemap.x_bins = ratemap.params.x;
ratemap.y_bins = ratemap.params.y;
ratemap.rateMaps = [];
save([basename '.ratemap.firingRateMap.mat'],'ratemap');

% 2- Detect (average) place fields
[placeFieldStats] = findPlaceFieldsAvg1D('firingMaps',firingMaps,'minPeak',1,'sepEdge',0.04);
save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');

% 3- Calculate template of place fields in the maze (for decoding, etc.)
[placeFieldTemplate] = findPlaceFieldsTemplate('placeFieldStats',placeFieldStats,'firingMaps',firingMaps);
save([basename '.placeFieldTemplate.mat'],'placeFieldTemplate');

% 4- Phase precession - NOT CHECKED
% theta phase
lfp = getLFP(refCh);
theta = bz_Filter(lfp,'passband',[5 15]);
% boundaries of each PF
for i=1:numel(spikes.UID) %
    for j=1:2
        for k=1:length(stats{i}{j}.peak)
            if stats{i}{j}.peak(k) ~= 0
                boundaries{i}{j}(k,1)= curve{i}{j}.x(stats{i}{j}.fieldX(k,1));
                boundaries{i}{j}(k,2)= curve{i}{j}.x(stats{i}{j}.fieldX(k,2));
            else
                boundaries{i}{j}(k,1)= NaN;
                boundaries{i}{j}(k,2)= NaN;
            end
        end
    end
end
% calculate phase precession
for i=1:numel(spikes.UID)
    for j=1:2
        for k=1:length(stats{i}{j}.x) %size(boundaries{i}{j},1)
            if ~isnan (stats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed
                count=count+1;
                [dataPP{count},statsPP{count}] = PhasePrecession(posTrials{j},spikes.times{i},theta.phase,'boundaries',boundaries{i}{j}(k,:));
                PlotPhasePrecession(dataPP{count},statsPP{count});
            end
        end
    end
end

