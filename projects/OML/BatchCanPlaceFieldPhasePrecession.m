function [points,smoothed,stats,basepath] = BatchCanPlaceFieldPhasePrecession(basepath)
% function [basepath] = BatchCanPlaceFieldPhasePrecession(basepath)

%%
if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

try
%     error
    load(fullfile(basepath,[basename '.phasePrecession.mat']),'PP');
    load(fullfile(basepath,[basename '.placeFields.cellinfo.mat']),'placeFieldStats');
    disp(['Session ' basepath ' already processed (phasePrecession file found).']);

    dataPP = PP.dataPP; statsPP = PP.statsPP;
    nBins = [10 25]*20; smooth = 0;
    for i=1:size(dataPP,1)
        for j=1:size(dataPP,2)
            if isempty(dataPP{i,j}), 
                points{i,j} = nan(0,2);
                continue; 
            end
            qq = [(dataPP{i,j}.position.x - statsPP{i,j}.boundaries(1))./diff(statsPP{i,j}.boundaries) dataPP{i,j}.position.phase];
            this = Restrict(qq,[0 1]); this(:,2) = wrap(this(:,2),2);
            points{i,j} = this;

            % Update stats because they were somehow computed taking into account spikes outside the field...?
            [beta,r,p] = CircularRegression(this(:,1),this(:,2),'slope',0);
            statsPP{i,j}.p(1) = p; statsPP{i,j}.slope(1) = beta(1); statsPP{i,j}.intercept(1) = r;
            
            this(:,1) = ceil(this(:,1)*nBins(1));
            this(:,2) = ceil(this(:,2)*nBins(2)/(2*pi)); this(this==0) = 1;
            saved{i,j} = Accumulate(this,1,'size',nBins)';
            smoothed{i,j} = Smooth(saved{i,j}/mean(saved{i,j}(:)),smooth,'type','cc');
            try
                stats{i,j} = [placeFieldStats.mapStats{i}{j}.x(1),placeFieldStats.mapStats{i}{j}.specificity(1), statsPP{i,j}.p(1), statsPP{i,j}.slope(1), statsPP{i,j}.intercept(1)];
            catch
                stats{i,j} = zeros(0,5);
            end
        end
    end
    

    %      clf
    %     x = linspace(0,1,nBins(1));
    %     y = linspace(0,4*pi,nBins(2)*2)/pi;
    %     for j=1:2
    %         q = cat(3,smoothed{:,j});
    %         subplot(1,2,j); PlotColorMap(repmat(nanmean(q,3),2,1),'x',x,'y',y);
    %     end
    %     title(basepath)
    %     drawnow

    return
end

%%

disp([datestr(clock) ': Starting session ' basepath]);

cd(basepath);

behavior = getStruct(basepath,'animal.behavior');
cell_metrics = getStruct(basepath,'cell_metrics');
spikesStruct = getStruct(basepath,'spikes.cellinfo');
try
    stim = behavior.stimON;
catch
   error(['No behavior.stimON found. Re-run script_OML_pipelinePlaceCells'])
end

try
    run = behavior.run;
catch
   error(['No behavior.run found. Re-run script_OML_pipelinePlaceCells'])
end

stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
nonstimIntervals = SubtractIntervals(run,stim); 

posTrials = behavior.positionTrialsRun;
if any(posTrials{1}(:,2)<0) || any(posTrials{1}(:,2)>1)
    error(['behavior.positionTrials are not normalized. Re-run script_OML_pipelinePlaceCells'])
end
d = cellfun(@(x) diff(x(:,2)),posTrials,'UniformOutput',false); d = cellfun(@(x) nanmedian(x(abs(x)>0)),d);
if any(d<0)
    error(['Some of behavior.positionTrials go backwards. Re-run script_OML_pipelinePlaceCells'])
end

%% Load data
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

%% 4- Phase precession
% theta phase

[parentFolder,dayName] = fileparts(basepath);
[~,animalName] = fileparts(parentFolder);
switch animalName
    case 'OML18'
        channel = 10;
    case 'OML19'
        channel = 34;
    case 'OLM21'
        channel = 22;
    case 'OML22'
        channel = 163;
    otherwise
        keyboard
end

lfp = getLFP(channel+1);
theta = bz_Filter(lfp,'passband',[5 15]);

% boundaries of each PF
for i=1:numel(spikesStruct.UID) %
    for j=1:2
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

phases = [theta.timestamps, theta.phase];
% calculate phase precession
tic
intervalsCell = {stimIntervals,nonstimIntervals};
for i=1:numel(spikesStruct.UID)
    for j=1:2
        intervals = intervalsCell{j};
        for k=1:length(placeFieldStats.mapStats{i}{j}.x) % number of place fields
            if ~isnan(placeFieldStats.mapStats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed
                this = posTrials{j}; 
                theseBoundaries = boundaries{i}{j}(k,:); theseBoundaries(theseBoundaries<0) = 0; theseBoundaries(theseBoundaries>1) = 1;
                if any(isnan(theseBoundaries)), continue; end
                if numel(Restrict(spikesStruct.times{i},intervals))<20, continue; end
                intervals = intervalsCell{j};
                [dataPP{i,j},statsPP{i,j}] = PhasePrecession(this,Restrict(spikesStruct.times{i},intervals),phases,'boundaries',theseBoundaries); 
            end
        end
    end
    disp(['Done with unit ' num2str(i) '/' num2str(numel(spikesStruct.UID)) '; estimated time remaining=' num2str(toc/i*(numel(spikesStruct.UID)-i))]);
end

PP = struct; PP.dataPP = dataPP; PP.statsPP = statsPP;
save(fullfile(basepath,[basename '.phasePrecession.mat']),'PP');

% for i=1:numel(dataPP)
%     if ~isempty(dataPP{i})
%      PlotPhasePrecession(dataPP{i},statsPP{i});
%      pause;
%      close all;
%     end
% end
% 
% nBins = [10 25];
% for i=1:size(dataPP,1)
%     for j=1:size(dataPP,2)
%         if isempty(dataPP{i,j}), continue; end
%         qq = [(dataPP{i,j}.position.x - statsPP{i,j}.boundaries(1))./diff(statsPP{i,j}.boundaries) dataPP{i,j}.position.phase];
%         this = Restrict(qq,[0 1]); this(:,2) = wrap(this(:,2),2);
%         this(:,1) = ceil(this(:,1)*nBins(1));
%         this(:,2) = ceil(this(:,2)*nBins(2)/(2*pi)); this(this==0) = 1;
%         saved{i,j} = Accumulate(this,1,'size',nBins)';
%         smoothed{i,j} = Smooth(saved{i,j}/sum(saved{i,j}(:)),1,'type','cc');
%     end
% end
% 
% figure(10); clf
% x = linspace(0,1,nBins(1));
% y = linspace(0,4*pi,nBins(2)*2)/pi;
% for j=1:2
%     q = cat(3,smoothed{:,j});
%     subplot(1,2,j); PlotColorMap(repmat(nanmean(q,3),2,1),'x',x,'y',y);
% end
% drawnow
% 
