function [dPeaks,dxCorr,dTheta,ratioCCG,specificity,dCenter,basepath] = BatchCanThetaCompression(basepath)

%%
if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,animalName] = fileparts(parentFolder);
sessionID = [animalName '_' dayName];

try
    load(fullfile(basepath,[basename '.thetaCompression.mat']),'thetaCompression');
    disp(['Session ' basepath ' already processed (thetaCompression file found).']);

    dPeaks = thetaCompression.dPeaks; dTheta = thetaCompression.dTheta; 
    ratioCCG = thetaCompression.ratioCCG; dxCorr = thetaCompression.dxCorr; 

    load(fullfile(basepath,[basename '.placeFields.cellinfo.mat']),'placeFieldStats');
    qq = cat(1,placeFieldStats.mapStats{:}); specificity = cellfun(@(x) x.specificity(1),qq);
    qq = cat(1,placeFieldStats.mapStats{:}); field = cellfun(@(x) [x.field false(100,1)],qq,'UniformOutput',0); field = cellfun(@(x) x(:,1),field,'UniformOutput',0);
    nUnits = size(placeFieldStats.mapStats,1); overlapping{1} = false(nUnits,nUnits); overlapping{2} = overlapping{1};
    dCenter = dPeaks;
    for i1 = 1:nUnits, for i2 = 1:nUnits, if i1~=i2,for j=1:2, overlapping{j}(i1,i2) = any(field{i1,j} & field{i2,j}); dCenter{j}(i1,i2) = mean(find(field{i1,j})) - mean(find(field{i2,j})); end; end; end; end
    return
end

disp([datestr(clock) ': Starting session ' basepath]);
%% Load place fields
behavior = getStruct(basepath,'animal.behavior');
run = behavior.run;
stim = behavior.stimON;

try
    cell_metrics = getStruct(basepath,'cell_metrics');
    spikesCell = cell_metrics.spikes.times';
catch
    error(['Error loading cell_metrics'])
end


try
    load(fullfile(basepath,[basename '.firingMapsAvg.cellinfo.mat']),'firingMaps');
catch
    error(['Firing maps not yet computed. Re-run script_OML_pipelinePlaceCells'])
end

try
    load(fullfile(basepath,[basename '.placeFields.cellinfo.mat']),'placeFieldStats');
catch
    error(['Place fields not yet computed. Re-run script_OML_pipelinePlaceCells'])
end

% if strcmp(animalName,'OML18')
%     cmPerTrack = nan(1,2);
% else
%     for j=1:2
%         from0to1 = behavior.positionTrials{j}; from0to1(isnan(from0to1(:,2)),:) = [];
%         from0to1 = from0to1(FindClosest(from0to1(:,2),[0.2 0.8]),:);
%         interpolated = interp1(behavior.timestamps(:),[behavior.position.x(:) behavior.position.y(:)],from0to1(:,1));
%         cmPerTrack(j) = sqrt(sum((interpolated(2,:)-interpolated(1,:)).^2,2))./diff(from0to1(:,2));
%     end
% end

nUnits = size(placeFieldStats.mapStats,1);

peakPosition = nan(nUnits,2);
for i=1:nUnits
    for j=1:2,
        peakPosition(i,j) = placeFieldStats.mapStats{i}{j}.x(1);
    end
end
nBins = length(firingMaps.rateMaps{1}{1});
cm = 300/nBins; 
peakPosition = peakPosition*cm; % convert from bins to cm (3m track).

dPeaks = cell(1,2); dxCorr = cell(1,2);
[dPeaks{1},dxCorr{1},dPeaks{2},dxCorr{2}] = deal(nan(nUnits,nUnits));
for i1=1:nUnits, for i2=1:nUnits, for j=1:2, if isnan(peakPosition(i1,j)) || isnan(peakPosition(i2,j)), continue; end
        dPeaks{j}(i1,i2) = peakPosition(i1,j)-peakPosition(i2,j);
        [x,xt] = xcorr(zscore(firingMaps.rateMaps{i1}{j}(:)),zscore(firingMaps.rateMaps{i2}{j}(:)),'coeff');
        xt = xt*cm;
        dxCorr{j}(i1,i2) = xt(findmax(x));
end; end; end

if size(behavior.trialID,2)>1, trialKind = 2-behavior.trialID(:,2); % trialKind=1 is for stim trials, trialKind=1 is for Off trials
else trialKind = behavior.trialID(:,1);
    if max(trialKind)==1, trialKind = 2-trialKind;
    end
end

spikesCellInField = cell(nUnits,2);
for i=1:nUnits
    for j=1:2
        if isnan(placeFieldStats.mapStats{i}{j}.x), continue; end
        in = InIntervals(spikesCell{i},SubtractIntervals(run,SubtractIntervals([0 Inf],behavior.trials(trialKind==j,:))));
        if j==1, in = in & InIntervals(spikesCell{i},stim); else, in = in & ~InIntervals(spikesCell{i},stim); end
        pos = interp1(behavior.positionTrials{j}(:,1),behavior.positionTrials{j}(:,2),spikesCell{i}); pos(~in) = nan;
        nBins = length(placeFieldStats.mapStats{i}{j}.field);
        pos = ceil(pos*nBins); pos(pos==0) = 1; pos(pos>nBins) = nBins;
        inField = false(size(pos));
        inField(~isnan(pos)) = placeFieldStats.mapStats{i}{j}.field(pos((~isnan(pos))));
        spikesCellInField{i,j} = spikesCell{i}(inField);
    end
end

dTheta = cell(1,2); [dTheta{1},dTheta{2}] = deal(nan(nUnits,nUnits));
[dThetaOpposite,ratioCCG,ratioCCGopposite] = deal(dTheta);
for i1=1:nUnits, for i2=1:nUnits, for j=1:2, if i1==i2 || isnan(peakPosition(i1,j)) || isnan(peakPosition(i2,j)), continue; end
            [h,ht] =  mPETH(spikesCellInField{i1,j},spikesCellInField{i2,j},'durations',[-1 1]*1,'nBins',2001);
            if ~any(h>0), continue; end
            f = FilterLFP([ht(:) h(:)],'passband',[2 30]);
            [h5,ht5] =  mPETH(spikesCellInField{i1,j},spikesCellInField{i2,j},'durations',[-1 1]*1,'nBins',401,'smooth',0); 
            h5 = h5.*length(spikesCellInField{i2,j}).*mode(diff(ht5(:))); % change from Hz to raw counts
            clear peaks; [peaks(:,1),peaks(:,2)] = SineWavePeaks(f,'mode','peaks'); peaks = Restrict(peaks,[-1 1]*0.1);
            if isempty(peaks), continue; end
            [~,index] = max(peaks(:,2)); peak = peaks(index,1);
            value = interp1(ht5(:),h5(:),peak(1)); value0 = min(h5(find(ht5>peak-0.075 & ht5<peak+0.075)));
            ok = value(1)>5;
            if ok, dTheta{j}(i1,i2) = peak(1); ratioCCG{j}(i1,i2) = value0/value; end
            % do the backup (opposite direction):
            peaks(sign(peaks(:,1))==sign(peak(1)),:) = []; % remove the peaks in the same direction, to locate the appropriate peak in the other direction
            if isempty(peaks), continue; end
            [~,index] = max(peaks(:,2)); peak = peaks(index,1);
            value = interp1(ht5(:),h5(:),peak(1)); value0 = min(h5(find(ht5>peak-0.075 & ht5<peak+0.075)));
            ok = value(1)>5;
            if ok, dThetaOpposite{j}(i1,i2) = peak(1); ratioCCGopposite{j}(i1,i2) = value0/value; end
end; end; end

for j=1:2
    l = lfit(dPeaks{j}(:),dTheta{j}(:));
    swap = abs(dTheta{j}-dPeaks{j}*l(1)+l(2)) > abs(dThetaOpposite{j}-dPeaks{j}*l(1)+l(2));
    dTheta{j}(swap) = dThetaOpposite{j}(swap); ratioCCG{j}(swap) = ratioCCGopposite{j}(swap);
end

thetaCompression = struct; 
thetaCompression.dPeaks = dPeaks; 
thetaCompression.dTheta = dTheta; 
thetaCompression.ratioCCG = ratioCCG;
thetaCompression.dxCorr = dxCorr; 
save(fullfile(basepath,[basename '.thetaCompression.mat']),'thetaCompression');
