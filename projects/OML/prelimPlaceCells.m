%% prelimnary place cell analysis Can's linear track experiment 

%%

% preprocessSession (label environment with gui_session)

% to get analog pulses, set manual threshold
[pulses] = getAnalogPulses('samplingRate',session.extracellular.sr,'manualThr',true);

% gather DLC tracking and sync with ttl; choose primary coordinates
general_behavior_file('primary_coords_dlc',4);

% coords = bestLED(behavior.position.x_1_point,behavior.position.y_1_point,behavior.position.x_2_point,behavior.position.y_2_point);

% get linear trial behavior (trials, running periods, etc); needs to set
% different speedTh for each session
linearTrackBehavior('manipulation',"true",'lapStart',15,'maze_sizes',300,'speedTh',20);

%%
% trials.outboundPos=Restrict([behavior.time' linpos],trials.outboundTs);
% trials.inboundPos=Restrict([behavior.time' linpos],trials.inboundTs);
% 
% posTrials{1}=trials.inboundPos;
% posTrials{2}=trials.outboundPos;

posTrials=behavior.positionTrialsRun;

%% Get theta LFP, phases, instantaneous frequency
thetaCh = 131;
lfp = getLFP(thetaCh); % 63 for day12; 24 for day13; 39 for day6; 46 for day7
%lfp.data = -lfp.data;% I'm taking the LFP from l-m and reverting its polarity
theta = bz_Filter(lfp,'passband',[5 15]);

%theta = Filter([lfp.timestamps lfp.data],'passband',[5 15]);
%hilb = hilbert(theta(:,2));
%phases = angle(hilb);
%phases = [lfp.timestamps, phases];
phases = [theta.timestamps, theta.phase];
lfpFreq= 1250/(2*pi)*diff(unwrap(phases)); 

save([basepath,filesep,[basename,'.theta.mat']],'phases','lfpFreq','thetaCh');
save([basepath,filesep,[basename,'.lfp.mat']],'lfp');
%% place cells
% ratemaps
spikes = importSpikes('CellType',"Pyramidal Cell"); % only load pyramidal neurons
[firingMaps] = firingMapAvg(posTrials,spikes,'nBins',100);

[placeFieldStats] = findPlaceFieldsAvg1D('firingMaps',firingMaps,'minPeak',1,'sepEdge',0.04);

[placeFieldTemplate] = findPlaceFieldsTemplate('placeFieldStats',placeFieldStats,'firingMaps',firingMaps);

%%
% place fields properties
for unit = 1:length(firingMaps.rateMaps)
    for t = 1:2
        if length(firingMaps.rateMaps{unit}) >= t && ~isempty(firingMaps.rateMaps{unit}{t}) %if having rate map for trial type t
           [info1,info2,spars,coef,select] = PlaceCellInfo(firingMaps.rateMaps{unit}{t},firingMaps.countMaps{unit}{t},firingMaps.occupancy{unit}{t});
           placeFieldStats.mapStats{unit}{t}.info1 = info1; clear info1;
           placeFieldStats.mapStats{unit}{t}.info2 = info2; clear info2;
           placeFieldStats.mapStats{unit}{t}.spars = spars; clear spars;
           placeFieldStats.mapStats{unit}{t}.coef = coef; clear coef;
           placeFieldStats.mapStats{unit}{t}.select = select; clear select;                      
        else
           placeFieldStats.mapStats{unit}{t}.info1 = NaN; 
           placeFieldStats.mapStats{unit}{t}.info2 = NaN; 
           placeFieldStats.mapStats{unit}{t}.spars = NaN;
           placeFieldStats.mapStats{unit}{t}.coef = NaN; 
           placeFieldStats.mapStats{unit}{t}.select = NaN; 
           placeFieldStats.mapStats{unit}{t}.peak(1) = NaN;
           placeFieldStats.mapStats{unit}{t}.mean(1) = NaN;
           placeFieldStats.mapStats{unit}{t}.size(1) = NaN;
        end
    end
end
for i = 1:2
    PCstats{i}.info1=[];PCstats{i}.info2=[];PCstats{i}.spars=[];PCstats{i}.specificity=[];PCstats{i}.coef=[];
    PCstats{i}.select=[];PCstats{i}.peak=[];PCstats{i}.mean=[];PCstats{i}.size=[];
end
for unit = 1:length(firingMaps.rateMaps)
    for i = 1:2
        if placeFieldStats.mapStats{unit}{i}.peak(1) > 0 % if peak firing rat larger than 0?
        PCstats{i}.info1 = cat(1,PCstats{i}.info1,placeFieldStats.mapStats{unit}{i}.info1);
        PCstats{i}.info2 = cat(1,PCstats{i}.info2,placeFieldStats.mapStats{unit}{i}.info2);
        PCstats{i}.spars = cat(1,PCstats{i}.spars,placeFieldStats.mapStats{unit}{i}.spars);
        PCstats{i}.specificity = cat(1,PCstats{i}.specificity,placeFieldStats.mapStats{unit}{i}.specificity);
        PCstats{i}.coef = cat(1,PCstats{i}.coef,placeFieldStats.mapStats{unit}{i}.coef);
        PCstats{i}.select = cat(1,PCstats{i}.select,placeFieldStats.mapStats{unit}{i}.select);
        PCstats{i}.peak = cat(1,PCstats{i}.peak,placeFieldStats.mapStats{unit}{i}.peak(1));
        PCstats{i}.mean = cat(1,PCstats{i}.mean,placeFieldStats.mapStats{unit}{i}.mean(1));
        PCstats{i}.size = cat(1,PCstats{i}.size,placeFieldStats.mapStats{unit}{i}.size(1));
        end
    end
end

save([basepath,filesep,[basename,'.PCstats.mat']],'PCstats');

%% Plotting place cell statistics
colorL = {'b','k'};
figure;
subplot(2,2,1);
for i =1:2
h = cdfplot(PCstats{i}.peak);hold on;set( h,'Color', colorL{i});
plot([nanmedian(PCstats{i}.peak) nanmedian(PCstats{i}.peak)],ylim,'-k');hold on;
end
p1=ranksum(PCstats{1}.peak,PCstats{2}.peak);
text(20,0.5,['p= ' num2str(p1,2)]);
title('peak FR');xlabel('firign rate (Hz)');ylabel('frac. of cells');xlim([-Inf 40]);
subplot(2,2,2);
for i =1:2
h=cdfplot(sqrt(PCstats{i}.size*2));set( h,'Color', colorL{i});hold on;
plot([nanmedian(sqrt(PCstats{i}.size*2)) nanmedian(sqrt(PCstats{i}.size*2))],ylim,'-k');hold on;
end
p1=ranksum(PCstats{1}.size*2,PCstats{2}.size*2);
text(5,0.5,['p= ' num2str(p1,2)]);
title('PF size');xlabel('cm');ylabel('frac. of cells');xlim([-Inf 15]);
subplot(2,2,3);
for i =1:2
h=cdfplot(PCstats{i}.info1);hold on;set( h,'Color', colorL{i});
plot([nanmedian(PCstats{i}.info1) nanmedian(PCstats{i}.info1)],ylim,'-k');hold on;
end
p1=ranksum(PCstats{1}.info1,PCstats{2}.info1);
text(1.5,0.5,['p= ' num2str(p1,2)]);
title('spat. info');xlabel('bits/spk');ylabel('frac. of cells');xlim([-Inf 3]);
subplot(2,2,4);
for i=1:2
h=cdfplot(PCstats{i}.spars);hold on;set( h,'Color', colorL{i});
plot([nanmedian(PCstats{i}.spars) nanmedian(PCstats{i}.spars)],ylim,'-k');hold on;
end
p1=ranksum(PCstats{1}.spars,PCstats{2}.spars);
text(0.6,0.5,['p= ' num2str(p1,2)]);
title('sparsity');xlabel('spars. coef.');ylabel('frac. of cells');
subtitle(['CA1 pyr n = ' num2str(length(PCstats{1}.peak)) '/ ' num2str(length(PCstats{2}.peak))]);

saveas(gcf,[basepath,filesep,'PCstats.png'],'png');

%% rate map for each lap - NOT FINISHED YET
load([basename '.animal.behavior.mat']);

firingMapsLaps = firingMapLaps(spikes,behavior,'nBins',100);

%save([basepath,filesep,[basename,'.firingMapsLaps.mat']],'firingMapsLaps');

%% unit osci freq - NEED CHECKING 

unitFreq = nan(length(spikes.UID),2); unitLFPfreq = nan(length(spikes.UID),2);
for u = 1:length(spikes.UID)
    for i = 1:2
        if placeFieldStats.mapStats{u}{i}.peak(1) ~= 0    
            spk = Restrict(spikes.times{u},trials{i}.timestampsRun); % restrict spike to only when running
            if i == 2 % i == 1; for stimulation trial number
            %spk = Restrict(spk,pulses.intsPeriods); % further restrict to when stimulation is On   
            spk = Restrict(spk,digitalIn.intsPeriods{7}); % for session with sinusoid
            end
            if length(spk) > 50 && strcmp(cell_metrics.putativeCellType{u},'Pyramidal Cell')
               tmp = Interpolate([phases(1:end-1,1) lfpFreq],spk);
               %lfpFreqSpk = mean(tmp(:,2));
               lfpFreqSpk = mean(tmp(:,3));
                [ccg2,ccg_time2] =CCG(spk,ones(length(spk)),'binSize',0.001,'duration',1);
                ccg_smth = Smooth(ccg2,10);
                
                [~,ind]=max(ccg_smth(551:651)); 
                if ind ~= 1 && ind ~= 101 
                   unitFreq(u,i) = 1/ccg_time2(551+ind);
                   unitLFPfreq(u,i) = unitFreq(u,i)-lfpFreqSpk;
                else
                   unitFreq(u,i)=NaN;  unitLFPfreq(u,i) = NaN;
                end
                
            else
                unitFreq(u,i)=NaN;    unitLFPfreq(u,i) = NaN;
            end
            clear spk ccg ind lfpFreqSpk tmp
        end
    end
end        

figure;
subplot(1,2,1);boxplot(unitFreq);title('unit freq');
subplot(1,2,2);boxplot(unitLFPfreq);title('unit freq - LFP freq');

saveas(gcf,[basepath,filesep,basename,'.unitFreq.png'],'png');

%% phase precession - NOT GOOD YET
 
% get boudaries from place fields that's normalized to 1
for i=1:numel(spikes.times) % for each pyramidal neuron
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

for i=1:numel(spikes.times)
    for j=1:2
        for k=1:length(placeFieldStats.mapStats{i}{j}.x) % number of place fields
            if ~isnan(placeFieldStats.mapStats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed
                %[dataPP{j}{i},statsPP{j}{i}] = PhasePrecession(posTrials{j},spikes.times{i},phases,'boundaries',boundaries{i}{j}(k,:));
                [dataPP{j}{i},statsPP{j}{i}] = fieldPhasePrecession(behavior.positionTrials{1,j},spikes.times{i},phases,'boundaries',boundaries{i}{j}(k,:));
            end
        end
    end
end

for i = 1:2
    PPstats{i}.slope=[]; PPstats{i}.intercept=[];PPstats{i}.r2=[]; PPstats{i}.p=[];
end
for i = 1:2
    for unit = 1:length(statsPP{i})
        if ~isempty(statsPP{i}{unit})
        PPstats{i}.slope = cat(1,PPstats{i}.slope,statsPP{i}{unit}.slope);
        PPstats{i}.intercept = cat(1,PPstats{i}.intercept,statsPP{i}{unit}.intercept);
        PPstats{i}.r2 = cat(1,PPstats{i}.r2,statsPP{i}{unit}.r2);
        PPstats{i}.p = cat(1,PPstats{i}.p,statsPP{i}{unit}.p);        
        end
    end
end

save([basepath,filesep,[basename,'.PP.mat']],'dataPP','statsPP','PPstats');

% plotting PP
for i=1:numel(dataPP)
    for j=1:numel(dataPP{i})
        if ~isempty(dataPP{i}{j})
            figure
         PlotPhasePrecession(dataPP{i}{j},statsPP{i}{j});
         pause;
         close all;
        end
    end
end

%% theta seq
basename = basenameFromBasepath(pwd);
load([basename '.animal.behavior.mat']);
load([basename '.ripples.events.mat']);

spk = importSpikes('CellType',"Pyramidal Cell");
spikesCell = spk.times';
spikes = Group(spikesCell{:});

% apply speed threshold to trial epochs
[quiet,quiescence] = QuietPeriods([behavior.timestamps' behavior.speed],0.1,0.5);
trials = SubtractIntervals(behavior.trials,quiet);

% or
trials = [trials{1,1}.timestampsRun;trials{1,2}.timestampsRun];
trials = sortrows(trials);

% Find theta cycles
try
    load([basepath,filesep,basename,'.thetaCyclesTask.mat'],'cycles');
catch
    lfpstruct = getLFP(thetaCh); % lfpstruct = lfp;
    lfp = lfpstruct.timestamps; lfp(:,2) = lfpstruct.data; % need to set it as timestamps first and load "data" later because "data" is an integer (and "lfp" should not be)
    [cycles,~] = FindThetaCycles(Restrict(lfp,trials));
    save([basepath,filesep,'.thetaCyclesTask.mat'],'cycles');
end
[windows] = SplitIntervals(cycles,'nPieces',6);
id = repmat((1:6)',length(cycles),1);

% calc theta seq

%     pos1 = behavior.positionTrialsRun{1};
%     t1 = behavior.positionTrialsRun{1}(:,1);
pos1 = [posTrials{1}(:,1) ZeroToOne(posTrials{1}(:,2))]; 
t1 = posTrials{1}(:,1);
intervals1 = t1(FindInterval(~isnan(pos1(:,2))));
in = repelem(InIntervals(cycles,stim),6,1);
[estimations1,actual1,errors1,average1] = ReconstructPosition(pos1,spikes,windows(~in,:),'training',intervals1,'id',id(~in));

%     pos2 = behavior.positionTrialsRun{2};
%     t2 = behavior.positionTrialsRun{2}(:,1);
pos2 = [posTrials{2}(:,1) ZeroToOne(posTrials{2}(:,2))];
t2 = posTrials{2}(:,1);
intervals2 = t2(FindInterval(~isnan(pos2(:,2))));
[estimations2,actual2,errors2,average2] = ReconstructPosition(pos2,spikes,windows(in,:),'training',intervals2,'id',id(in));

% Plotting
figure
subplot(1,2,1);
PlotColorMap([average1 average1]);
PlotHVLines(100,'h','w--','linewidth',2);
[~,~,a,b] = FindReplayScore(average1,'circular','off');
plot([1 6],[a b],'k');
plot([1 6]+6,[a b],'k');ylim([60 140]);
title([num2str(a) ' ' num2str(b)]);
set(gca,'xtick',[]);
xlabel('theta phase');title('stim OFF');

subplot(1,2,2);
PlotColorMap([average2 average2]);
PlotHVLines(100,'h','w--','linewidth',2);
[~,~,a,b] = FindReplayScore(average2,'circular','off');
plot([1 6],[a b],'k');
plot([1 6]+6,[a b],'k');ylim([60 140]);
title([num2str(a) ' ' num2str(b)]);
set(gca,'xtick',[]);
xlabel('theta phase');title('stim ON');                

%
% saveas(gcf,[basepath,filesep,basename,'.thetaSequenceDecoding.png'],'png');    
%% Plotting raster plot for all laps for each individual neuron

% Loading
basepath = pwd;
basename = basenameFromBasepath(basepath);
load([basepath,filesep,[basename,'.spikes.cellinfo.mat']]);
load([basepath,filesep,[basename,'.animal.behavior.mat']]);
spikes = importSpikes('CellType',"Pyramidal Cell"); % only load pyramidal neurons

%
for i=1:numel(spikes.times) % for each py neuron
    for c=1  % for each condition

        unitFiring = Restrict(spikes.times{1,i},trials{1,c}.timestampsRun);
        positions = behavior.positionTrialsRun{1,c}(:,2);
        posTimestamps = behavior.positionTrialsRun{1,c}(:,1);

        pos = interp1(posTimestamps,positions,unitFiring); % interpolate animal position while cell fired
        speed = interp1(posTimestamps,positions,unitFiring); 
        [~,trialID] = InIntervals(unitFiring,trials{1,c}.timestamps); % find trialID of each spike
        ok = trialID>0; % trialID = 0 means firing was not in any trial

        scatter(pos,unitFiring);hold on
        plot(positions,posTimestamps)
        
        %{
        for i=1:size(trials,1)
        curve(i,:) = FiringCurve(Restrict([posTimestamps, positions],trials(i,:)),Restrict(unitFiring,trials(i,:)));
        end
        %}

        % extract rmap for each lap from structure
        rmapLap = [];
        for t = 1:size(firingMapsLaps.rateMaps{1,1},1)
            rmapLap = [rmapLap;firingMapsLaps.rateMaps{i,1}{t,c}];
        end

        subplot(1,2,1);
        PlotColorMap(rmapLap);
        subplot(1,2,2);
        RasterPlot([pos(ok) trialID(ok)]);
        ylim([0 max(trialID)+2])% generate Raly style raster plot

    end
end

