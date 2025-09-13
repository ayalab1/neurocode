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

%% Loading posTrials & spikes
load('posTrials.mat')
basename = basenameFromBasepath(pwd);
load(fullfile(pwd,[basename,'.spikes.cellinfo.mat']))

%% 1-Calculate firing maps
spikes = importSpikes('CellType',"Pyramidal Cell");
[firingMaps] = firingMapAvg(posTrials,spikes);
% save([basename '.firingMapsAvg.cellinfo.mat'],'firingMaps');

%% 1.1- Save struct for cell explorer (you can display ratemaps in cell explorer)
% Optimized for only 1d, though 2d works too
ratemap = firingMaps;
for i = 1:length(ratemap.rateMaps)
    ratemap.map{i} = ratemap.rateMaps{i}{1};
end
ratemap.x_bins = ratemap.params.x;
ratemap.y_bins = ratemap.params.y;
ratemap.rateMaps = [];
save([basename '.ratemap.firingRateMap.mat'],'ratemap');

%% 2- Detect (average) place fields
[placeFieldStats] = findPlaceFieldsAvg1D('firingMaps',firingMaps,'minPeak',1,'sepEdge',0.04);
% save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');

% 2.1 alternative to detect place fields 2D
for i = 1:length(firingMaps.rateMaps)
    map.z = firingMaps.rateMaps{i}{1};
    map.time = firingMaps.occupancy{i}{1};
    map.count = firingMaps.countMaps{i}{1};
    map.x = firingMaps.params.x;
    map.y = firingMaps.params.y;
    placeFieldStats{i} = MapStats(map);
end
save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');

%% 3- Calculate template of place fields in the maze (for decoding, etc.)
[placeFieldTemplate] = findPlaceFieldsTemplate('placeFieldStats',placeFieldStats,'firingMaps',firingMaps);
% save([basename '.placeFieldTemplate.mat'],'placeFieldTemplate');

%% 4- Phase precession - NOT CHECKED
% theta phase
lfp = getLFP(refCh);
theta = bz_Filter(lfp,'passband',[5 15]);

%% boundaries of each PF
for i=1:numel(spikes.UID) %
    for j=1%:2
        for k=1:length(placeFieldStats.mapStats{i}{j}.peak)
            if placeFieldStats.mapStats{i}{j}.peak(k) ~= 0
                boundaries{i}{j}(k,1)= firingMaps.params.x(placeFieldStats.mapStats{i}{j}.fieldX(k,1));
                boundaries{i}{j}(k,2)= firingMaps.params.x(placeFieldStats.mapStats{i}{j}.fieldX(k,2));
            else
                boundaries{i}{j}(k,1)= NaN;
                boundaries{i}{j}(k,2)= NaN;
            end
        end
    end
end

count = 0;
phases = [theta.timestamps, theta.phase];
% calculate phase precession
for i=1:numel(spikes.UID)
    for j=1%:2
        for k=1:length(placeFieldStats.mapStats{i}{j}.x) % number of place fields
            if ~isnan(placeFieldStats.mapStats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed
                count=count+1;
                [dataPP{count},statsPP{count}] = PhasePrecession(posTrials{j},spikes.times{i},phases,'boundaries',boundaries{i}{j}(k,:));
            end
        end
    end
end

PP = struct; PP.dataPP = dataPP; PP.statsPP = statsPP;
%save([basename '.phasePrecession.mat'],'PP');

for i=1:numel(dataPP)
    if ~isempty(dataPP{i})
        figure
     PlotPhasePrecession(dataPP{i},statsPP{i});
     %pause;
     %close all;
    end
end

