%% Multi-PSTH
% Plot multiple PSTHs for multiple regions, barrage vs ripples

function BARR_PSTH(plotSave,useState)
% plotSave should be char of where you want the plots saved
% useState should be string of the state you want to look at
if nargin <2
    useState = "NREM";
    warning('No state specified, defaulting to NREM');
end
close all
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
HSE = loadBarrage(basepath);
if ~isempty(HSE.NREM)
load(strcat(basename,'.cell_metrics.cellinfo.mat'));
load(strcat(basepath,filesep,basename,'.SleepState.states.mat'));
switch useState
    case "nonTheta"
        if ~exist([basepath filesep 'Barrage_Files' filesep 'nonTheta'])
            mkdir([basepath filesep 'Barrage_Files' filesep 'nonTheta']);
        end
        if contains(plotSave,'clean')
            plotSave = [basepath filesep 'Barrage_Files' filesep 'nonTheta' filesep basename '.clean.'];
        else
            plotSave = [basepath filesep 'Barrage_Files' filesep 'nonTheta' filesep basename '.'];
        end
    case "imm"
        if ~exist([basepath filesep 'Barrage_Files' filesep 'imm'])
            mkdir([basepath filesep 'Barrage_Files' filesep 'imm']);
        end
        if contains(plotSave,'clean')
            plotSave = [basepath filesep 'Barrage_Files' filesep 'imm' filesep basename '.clean.'];
        else
            plotSave = [basepath filesep 'Barrage_Files' filesep 'imm' filesep basename '.'];
        end
    case "theta"
        if ~exist([basepath filesep 'Barrage_Files' filesep 'theta'])
            mkdir([basepath filesep 'Barrage_Files' filesep 'theta']);
        end
        if contains(plotSave,'clean')
            plotSave = [basepath filesep 'Barrage_Files' filesep 'theta' filesep basename '.clean.'];
        else
            plotSave = [basepath filesep 'Barrage_Files' filesep 'theta' filesep basename '.'];
        end
    case "REM"
        if ~exist([basepath filesep 'Barrage_Files' filesep 'REM'])
            mkdir([basepath filesep 'Barrage_Files' filesep 'REM']);
        end
        if contains(plotSave,'clean')
            plotSave = [basepath filesep 'Barrage_Files' filesep 'REM' filesep basename '.clean.'];
        else
            plotSave = [basepath filesep 'Barrage_Files' filesep 'REM' filesep basename '.'];
        end
    otherwise
        useState = "NREM";
end

if ~exist([basepath filesep 'Barrage_Files' filesep 'PSTHmet'])
    mkdir([basepath filesep 'Barrage_Files' filesep 'PSTHmet']);
end
if ~exist(strcat(basepath,filesep,basename,'.ripples.events.mat'))
    if exist(strcat(basepath,filesep,basename,'SWR.events.mat'))
        load(strcat(basepath,filesep,basename,'SWR.events.mat'));
        ripples = SWR;
        if isfield(ripples,'peaktimes')
            ripples.peaks = ripples.peaktimes;
        end
        plotRips = 1;
    elseif exist(strcat(basepath, filesep, 'old files', filesep, 'SWR.mat'))
        load(strcat(basepath, filesep, 'old files', filesep, 'SWR.mat'));
        ripples = SWR;
        if isfield(ripples,'peaktimes')
            ripples.peaks = ripples.peaktimes;
        end
        plotRips = 1;
    elseif exist(strcat(basepath, filesep, 'old files', filesep, 'SWR.mat'))
        load(strcat(basepath, filesep, 'old files', filesep, 'SWR.mat'));
        ripples = SWR;
        if isfield(ripples,'peaktimes')
            ripples.peaks = ripples.peaktimes;
        end
        plotRips = 1;
    else
        plotRips = 0;
    end
else
    load(strcat(basename,'.ripples.events.mat'));
    if exist('SWR')
        ripples = SWR;
    end
    plotRips = 1;
end
if plotRips
    if isfield(ripples, 'peaktimes')
        ripples.peaks = ripples.peaktimes;
    end
end
%%

if useState == "nonTheta"
    [HSEuse.peaks,useIdx] = Restrict(HSE.peaks,SleepState.ints.nonTHETA);
    HSEuse.timestamps = HSE.timestamps(useIdx,:);
    HSEuse.amplitudes = HSE.amplitudes(useIdx);
    HSEuse.duration = HSE.duration(useIdx);
elseif useState == "imm"
    HSEuse.peaks = HSE.peaks(HSE.keep(HSE.imm));
    HSEuse.timestamps = HSE.timestamps(HSE.keep(HSE.imm),:);
    HSEuse.amplitudes = HSE.amplitudes(HSE.keep(HSE.imm));
    HSEuse.duration = HSE.duration(HSE.keep(HSE.imm));
elseif useState == "theta"
    [wake_theta] = getWakeStates(SleepState);
    if ~isempty(wake_theta)
        [HSEuse.peaks,useIdx] = Restrict(HSE.peaks,wake_theta);
        HSEuse.timestamps = HSE.timestamps(useIdx,:);
        HSEuse.amplitudes = HSE.amplitudes(useIdx);
        HSEuse.duration = HSE.duration(useIdx);
    else
        warning('no theta, returning');
        return
    end
elseif useState == "REM"
    [HSEuse.peaks,useIdx] = Restrict(HSE.peaks,SleepState.ints.REMstate);
    HSEuse.timestamps = HSE.timestamps(useIdx,:);
    HSEuse.amplitudes = HSE.amplitudes(useIdx);
    HSEuse.duration = HSE.duration(useIdx);
else
    if useState ~= "NREM"
        warning('Unrecognized state, defaulting to NREM');
    end
    if isfield(HSE, 'keep')
        HSEuse.timestamps = HSE.timestamps(HSE.keep,:);
        HSEuse.peaks = HSE.peaks(HSE.keep);
%         HSEuse.amplitudes = HSE.amplitudes(HSE.keep);
        HSEuse.duration = HSE.duration(HSE.keep);
    else
        HSEuse = eventIntervals(HSE,SleepState.ints.NREMstate,1);
    end
    if isfield(HSE,'keep')&&isfield(HSE,'NREM')
        HSEuse.timestamps = HSEuse.timestamps(HSE.NREM,:);
        HSEuse.peaks = HSEuse.peaks(HSE.NREM);
%         HSEuse.amplitudes = HSEuse.amplitudes(HSE.NREM);
        HSEuse.duration = HSEuse.duration(HSE.NREM);
    else
        HSEuse = eventIntervals(HSEuse,SleepState.ints.NREMstate,1);
    end
end

if isempty(HSEuse.timestamps)
    warning('no events during this sleep state, returning');
    return
end


%% PSTH
tempUseReg = unique(cell_metrics.brainRegion);
regions = []; tuC = 1; y1 = 0; y2 = 0; y3=0;
for i = 1:length(tempUseReg)
    if contains(convertCharsToStrings(tempUseReg{i}),"CA1")
        if ~y1
            regions{tuC} = "CA1";
            tuC = tuC+1; y1 = 1;
        end
    elseif contains(convertCharsToStrings(tempUseReg{i}),"CA2")
        if (~y2)
            regions{tuC} = "CA2";
            tuC = tuC+1; y2 = 1;
        end
    elseif contains(convertCharsToStrings(tempUseReg{i}),"CA3")
        if (~y3)
            regions{tuC} = "CA3";
            tuC = tuC+1; y3 = 1;
        end
    else
        regions{tuC} = convertCharsToStrings(tempUseReg{i});
        tuC = tuC+1;
    end
end

check = ["CA1" "CA2" "CA3" "CTX" "DG" "MEC" "LEC"];
regTot = 0;
for i = 1:length(check)
    for j = 1:length(regions)
        regCheck = regions{j};
        if contains(regCheck,check(i))
            regTot = regTot+1;
        end
    end
end

if (~isfield(cell_metrics.tags,'Pb'))&&(~isfield(cell_metrics.tags,'Nb'))&&(~isfield(cell_metrics.tags,'Pr'))&&(~isfield(cell_metrics.tags,'Nr'))
    runTags = 1;
else
    runTags = 0;
end
if runTags
    if ~isfield(cell_metrics.tags,'Pb')
        cell_metrics.tags.Pb = [];
    end
    if ~isfield(cell_metrics.tags,'Nb')
        cell_metrics.tags.Nb = [];
    end
    if ~isfield(cell_metrics.tags,'Pr')
        cell_metrics.tags.Pr = [];
    end
    if ~isfield(cell_metrics.tags,'Nr')
        cell_metrics.tags.Nr = [];
    end
end

subNum = 1;
for i = 1:length(check)
    for j=1:length(regions)
        regCheck = regions{j};
        if contains(regCheck, check(i))
            br = convertStringsToChars(check(i));
            if ~exist([basepath filesep 'Barrage_Files' filesep basename '.' br 'pyr.cellinfo.mat'])
                spikes = importSpikes('brainRegion', check(i), 'cellType', "Pyramidal Cell");
                save([basepath filesep 'Barrage_Files' filesep basename '.' br 'pyr.cellinfo.mat'], 'spikes');
            else
                load([basepath filesep 'Barrage_Files' filesep basename '.' br 'pyr.cellinfo.mat']);
            end
            if ~isempty(spikes.times)
                if plotRips
                    %                     PSTH_ripples = computePSTH(ripples,spikes,'duration',1,'plots', false,'alignment','peaks','zscorePlot',false);
%                     [PSTH_ripples,PSTH_t] = makePETH(spikes,ripples.timestamps(:,1),[-3 3]);
                    [PSTH_ripples, PSTH_t] = manualPSTH(spikes, ripples.timestamps(:,1),[-1,1],1/1000);
                    if runTags
                        warning('NOT ADJUSTED TO NEW DIMENSIONS');
%                         for k = 1:size(PSTH_ripples,1)
%                            if ((max(PSTH_ripples(k,91:111))-max(PSTH_ripples(k,1:71))) > 1.5*std(PSTH_ripples(k,:)))
%                                cell_metrics.tags.Pr = [cell_metrics.tags.Pr spikes.UID(k)];
%                            else
%                                cell_metrics.tags.Nr = [cell_metrics.tags.Nr spikes.UID(k)];
%                            end
%                         end
                    end
%                     [PSTH_ripples,PSTH_t] = makePETH(spikes,ripples.peaks,[-1 1]);
%                     index = getSort(PSTH_ripples);
                    index = 1:size(PSTH_ripples,1);
%                     PSTH_ripples = manualZscore(PSTH_ripples,spikes);
%                     PSTH_ripples = zscore(PSTH_ripples,1,2);
                    plotPSTH(PSTH_ripples, PSTH_t, index, regTot, subNum, strcat(check(i), '/ripples'),1);
                    PSTHmets.PSTH_ripples = PSTH_ripples;
                    PSTHmets.PSTH_ripT = PSTH_t;
                    PSTHmets.UID = spikes.UID;
                end
                %                     PSTH_bar = computePSTH(HSEnrem,spikes,'duration',1,'plots', false,'alignment','peaks','zscorePlot',false);
%                 [PSTH_bar,PSTH_t] = makePETH(spikes,HSEuse.timestamps(:,1),[-3 3]);
                [PSTH_bar, PSTH_t] = manualPSTH(spikes, HSEuse.timestamps(:,1),[-1 1],1/1000);
                if runTags
                    warning('NOT ADJUSTED TO NEW DIMENSIONS');
%                     for k = 1:size(PSTH_bar,1)
%                         if (mean(PSTH_bar(k,1001:2001))-mean(PSTH_bar(k,1:100))) > 0.3
%                             cell_metrics.tags.Pb = [cell_metrics.tags.Pb spikes.UID(k)];
%                         else
%                             cell_metrics.tags.Nb = [cell_metrics.tags.Nb spikes.UID(k)];
%                         end
%                     end
                end
%                 PSTH_bar = manualZscore(PSTH_bar,spikes);
%                 PSTH_bar = zscore(PSTH_bar,1,2);
%                 [PSTH_bar,PSTH_t] = makePETH(spikes,HSEuse.peaks,[-1 1]);
                plotPSTH(PSTH_bar, PSTH_t, index, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'),1);
                PSTHmets.PSTH_bar = PSTH_bar;
                PSTHmets.PSTH_barT = PSTH_t;

                if useState == "NREM"
                    if contains(plotSave, 'clean')
                        save([basepath filesep 'Barrage_Files' filesep 'PSTHmet' filesep basename '.clean.' br 'PSTHmets.mat'], 'PSTHmets');
                    else
                        save([basepath filesep 'Barrage_Files' filesep 'PSTHmet' filesep basename '.' br 'PSTHmets.mat'], 'PSTHmets');
                    end
                else
                    if contains(plotSave, 'clean')
                        save(strcat(basepath, filesep, 'Barrage_Files', filesep, useState, filesep, basename, '.clean.', br, 'PSTHmets.mat'), 'PSTHmets');
                    else
                        save(strcat(basepath, filesep, 'Barrage_Files', filesep, useState, filesep, basename, '.', br, 'PSTHmets.mat'), 'PSTHmets');
                    end
                end

            else
                warning(['No pyramidal cells in ' br]);
            end

            subNum = subNum+1; % should max at 4
        end
    end
end
sgtitle('Pyramidal PSTH');
saveas(gcf,[plotSave 'PSTH.png']);

subNum = 1; intEmpt = 0;
for i = 1:length(check)
    for j=1:length(regions)
        regCheck = regions{j};
        if contains(regCheck, check(i))
            br = convertStringsToChars(check(i));
            spikes = [];
            if ~exist([basepath filesep 'Barrage_Files' filesep basename '.' br 'int.cellinfo.mat'])
                spikes = importSpikes('brainRegion', check(i), 'cellType', ["Wide Interneuron"; "Narrow Interneuron"]);
                save([basepath filesep 'Barrage_Files' filesep basename '.' br 'int.cellinfo.mat'], 'spikes');
            else
                load([basepath filesep 'Barrage_Files' filesep basename '.' br 'int.cellinfo.mat']);
            end
            if ~isempty(spikes.times)
                if plotRips
                    %                     PSTH_ripples = computePSTH(ripples,spikes,'duration',1,'plots', false,'alignment','peaks','zscorePlot',false);
%                     [PSTH_ripples,PSTH_t] = makePETH(spikes,ripples.timestamps(:,1),[-3 3]);
%                     [PSTH_ripples,PSTH_t] = makePETH(spikes,ripples.peaks,[-1 1]);
                    [PSTH_ripples, PSTH_t] = manualPSTH(spikes, ripples.timestamps(:,1),[-1 1],1/1000);
                    if runTags
                        warning('NOT ADJUSTED TO NEW DIMENSIONS');
%                         for k = 1:size(PSTH_ripples,1)
%                             if (max(PSTH_ripples(k,91:111))-max(PSTH_ripples(k,1:71))) > 1.5*std(PSTH_ripples(k,:))
%                                 cell_metrics.tags.Pr = [cell_metrics.tags.Pr spikes.UID(k)];
%                             else
%                                 cell_metrics.tags.Nr = [cell_metrics.tags.Nr spikes.UID(k)];
%                             end
%                         end
                    end
%                     index = getSort(PSTH_ripples);
                    index = 1:size(PSTH_ripples,1);
%                     PSTH_ripples = manualZscore(PSTH_ripples,spikes);
%                     PSTH_ripples = zscore(PSTH_ripples,1,2);
                    plotPSTH(PSTH_ripples,  PSTH_t, index, regTot, subNum, strcat(check(i), '/ripples'),2);
                    PSTHmetsINT.PSTH_ripples = PSTH_ripples;
                    PSTHmetsINT.ripT = PSTH_t;
                    PSTHmetsINT.UID = spikes.UID;
                end
                %                 PSTH_bar = computePSTH(HSEnrem,spikes,'duration',1,'plots', false,'alignment','peaks','zscorePlot',false);
%                 [PSTH_bar,PSTH_t] = makePETH(spikes,HSEuse.timestamps(:,1),[-3 3]);
                [PSTH_bar, PSTH_t] = manualPSTH(spikes, HSEuse.timestamps(:,1),[-1 1],1/1000);
                if runTags
                    warning('NOT ADJUSTED TO NEW DIMENSIONS');
%                     for k = 1:size(PSTH_bar,1)
%                         if (mean(PSTH_bar(k,101:201))-mean(PSTH_bar(k,1:100))) > 0.3
%                             cell_metrics.tags.Pb = [cell_metrics.tags.Pb spikes.UID(k)];
%                         else
%                             cell_metrics.tags.Nb = [cell_metrics.tags.Nb spikes.UID(k)];
%                         end
%                     end
                end
%                 [PSTH_bar,PSTH_t] = makePETH(spikes,HSEuse.peaks,[-1 1]);
%                 PSTH_bar = manualZscore(PSTH_bar,spikes);
%                 PSTH_bar = zscore(PSTH_bar, 1, 2);
                plotPSTH(PSTH_bar, PSTH_t, index, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'),2);
                PSTHmetsINT.PSTH_bar = PSTH_bar;
                PSTHmetsINT.barT = PSTH_t;

                if useState == "NREM"
                    if contains(plotSave, 'clean')
                        save([basepath filesep 'Barrage_Files' filesep 'PSTHmet' filesep basename '.clean.' br 'PSTHmetsINT.mat'], 'PSTHmetsINT');
                    else
                        save([basepath filesep 'Barrage_Files' filesep 'PSTHmet' filesep basename '.' br 'PSTHmetsINT.mat'], 'PSTHmetsINT');
                    end
                else
                    if contains(plotSave, 'clean')
                        save(strcat(basepath, filesep, 'Barrage_Files', filesep, useState, filesep, basename, '.clean.', br, 'PSTHmetsINT.mat'), 'PSTHmetsINT');
                    else
                        save(strcat(basepath, filesep, 'Barrage_Files', filesep, useState, filesep, basename, '.', br, 'PSTHmetsINT.mat'), 'PSTHmetsINT');
                    end
                end

            else
                warning(['No interneurons in ' br]);
                intEmpt = intEmpt+1;
            end
            subNum = subNum+1; % should max at 4
        end
    end
end
if ~(intEmpt == length(regions))
    sgtitle('Interneuron PSTH');
    saveas(gcf,[plotSave 'PSTHint.png']);
end

if runTags
    warning('NOT ADJUSTED TO NEW DIMENSIONS');
%     cell_metrics.tags.Pb = unique(cell_metrics.tags.Pb);
%     cell_metrics.tags.Nb = unique(cell_metrics.tags.Nb);
%     cell_metrics.tags.Pr = unique(cell_metrics.tags.Pr);
%     cell_metrics.tags.Nr = unique(cell_metrics.tags.Nr);
%     save(strcat(basename,'.cell_metrics.cellinfo.mat'), 'cell_metrics');
end


%% Event duration
duration = HSEuse.timestamps(:,2)-HSEuse.timestamps(:,1);
figure(3);
hist(duration,50);title('barrages');xlabel('duration (s)');
saveas(gcf,[plotSave 'evtDist.png']);

%% CCG
if plotRips
%         t_ripple{1}=ripples.peaks(:);
%         t_barrages{1}=HSEnrem.peaks(:);
    t_ripple{1} = ripples.timestamps(:,1);
    t_barrages{1} = HSEuse.timestamps(:,1);
    t_ripple_id{1}=ones(1,length(ripples.peaks(:)))';
    t_barrage_id{1}=ones(1,length(HSEuse.peaks(:)))';
    binsize=0.1;duration=20;
    %ccg(t,i,j) number of events of group j at time lag t with respect to
    %reference events from group i
    [ccg_ripple_barrage,t_ripple_barrage] = CCG(cat(1,t_ripple{1},t_barrages{1}),cat(1,t_ripple_id{1},2*t_barrage_id{1}),'binSize',binsize,'duration',duration,'norm','rate');
        figure(4);plot(t_ripple_barrage,ccg_ripple_barrage(:,2,1),'r');hold on;title('ccg barr-SWR'); yline(0,'k--');xlim([-3 3]);
        saveas(gcf,[plotSave 'CCG.png']);
% Raly shuffle
r = ripples.timestamps(:,1);
iri = diff([0;r]);
[~,randomOrder] = sort(rand(size(iri)));
shuffled_iri = iri(randomOrder);
shuffled_r = cumsum([shuffled_iri]); %raly had cumsum([0; shuffled_iri]);

[ccg_ripple_barrage_shuf,t_ripple_barrage_shuf] = CCG(cat(1,shuffled_r,t_barrages{1}),cat(1,t_ripple_id{1},2*t_barrage_id{1}),'binSize',binsize,'duration',duration,'norm','rate');
        figure(5);plot(t_ripple_barrage_shuf,ccg_ripple_barrage_shuf(:,2,1),'r');hold on;title('ccg barr-SWR'); yline(0,'k--');xlim([-3 3]);
        saveas(gcf,[plotSave 'CCG_shuff.png']);
        % Old jitter start
%     for it = 1:10000
%         t_jitter = HSEuse.timestamps(:,1);
%         for j = 1:100
%             t_jitter = t_jitter + (5-(-5))*rand()-5;
%         end
%         [temp_ripple_jitter,t_ripple_jitter] = CCG(cat(1,t_ripple{1},t_jitter),cat(1,t_ripple_id{1},2*t_barrage_id{1}),'binSize',binsize,'duration',duration,'norm','rate');
%         if it==1
%             ccg_ripple_jitter = temp_ripple_jitter(:,2,1);
%             ccg_ripple_jitter_rev = temp_ripple_jitter(:,1,2);
%         else
%             ccg_ripple_jitter(:,it) = temp_ripple_jitter(:,2,1);
%             ccg_ripple_jitter_rev(:,it) = temp_ripple_jitter(:,1,2);
%         end
%     end
    % Old jitter end
    %     figure(5);plot(t_ripple_jitter,mean(ccg_ripple_jitter,2),'r');hold on;title('ccg jitter-SWR');
    %     saveas(gcf,[plotSave 'CCGjitter.png']);
    
    if useState == "NREM"
        CCG_dat.time = t_ripple_barrage; CCG_dat.y = ccg_ripple_barrage(:,2,1); CCG_dat.shuffle = ccg_ripple_barrage_shuf(:,2,1); CCG_dat.jitterT=t_ripple_barrage_shuf;
        CCG_dat.reverse = ccg_ripple_barrage(:,1,2); CCG_dat.shuffle_rev = ccg_ripple_barrage_shuf(:,1,2);
        save([plotSave 'CCG_dat.mat'],'CCG_dat');
    elseif useState == "imm"
        CCG_dat.time = t_ripple_barrage; CCG_dat.y = ccg_ripple_barrage(:,2,1); CCG_dat.shuffle = ccg_ripple_barrage_shuf(:,2,1); CCG_dat.jitterT=t_ripple_barrage_shuf;
        CCG_dat.reverse = ccg_ripple_barrage(:,1,2); CCG_dat.shuffle_rev = ccg_ripple_barrage_shuf(:,1,2);
        save([plotSave 'CCG_dat.mat'],'CCG_dat');
    end
end

%% Fraction per state
if useState == "NREM"
    [wake_theta] = getWakeStates(SleepState);
    num_theta = 0; num_nontheta = 0; num_REM = 0; num_NREM = 0;
    if ~isempty(wake_theta)
        temp_theta = Restrict(HSE.peaks(HSE.keep),wake_theta);
    else
        temp_theta = [];
    end
    num_theta = length(temp_theta);
    temp_wake = Restrict(HSE.peaks(HSE.keep),SleepState.ints.WAKEstate);
    num_nontheta = length(temp_wake) - num_theta; %if it's wake but not theta, must be nontheta wake
    num_REM = length(Restrict(HSE.peaks(HSE.keep),SleepState.ints.REMstate));
    num_NREM = length(Restrict(HSE.peaks(HSE.keep),SleepState.ints.NREMstate));
    if isfield(SleepState.ints,'THETA')
        timeTheta =sum(SleepState.ints.THETA(:,2)-SleepState.ints.THETA(:,1));
    else
        timeTheta = 1;
    end
    if isfield(SleepState.ints,'nonTHETA')
        timeNonTheta = sum(SleepState.ints.nonTHETA(:,2)-SleepState.ints.nonTHETA(:,1));
    else
        timeNonTheta = 1;
    end
    if isfield(SleepState.ints,'REMstate')
        timeREM = sum(SleepState.ints.REMstate(:,2)-SleepState.ints.REMstate(:,1));
    else
        timeREM = 1;
    end
    if isfield(SleepState.ints,'NREMstate')
        timeNREM = sum(SleepState.ints.NREMstate(:,2)-SleepState.ints.NREMstate(:,1));
    else
        timeNREM = 1;
    end
    figure(6); bar([num_theta, num_nontheta, num_REM, num_NREM]);
    xticklabels({'Theta', 'Nontheta', 'REM', 'NREM'});
    title('Distribution of events across states');
    ylabel('Number of events');
    labelsState = string([num_theta, num_nontheta, num_REM, num_NREM]);
    text([1 2 3 4],[num_theta, num_nontheta, num_REM, num_NREM],labelsState,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    saveas(gcf,[plotSave 'stateSplit.png']);
    stateCount.theta = num_theta; stateCount.nonTheta = num_nontheta; stateCount.REM = num_REM; stateCount.NREM = num_NREM;
    stateCount.timeTheta = timeTheta; stateCount.timeNonTheta = timeNonTheta; stateCount.timeREM = timeREM; stateCount.timeNREM = timeNREM;
    save([plotSave 'stateCount.mat'],'stateCount');
end
%% Fraction per task state
% maybe we could change this to be a stacked bar plot per task portion (so
% pre on x, then stack wake/nontheta/rem/nrem
if useState == "NREM" %only need this for main
    load([basename '.session.mat']);
    sessEp = cell(1, length(session.epochs));
    xtickuse = cell(1, length(session.epochs));
    bary = nan(1, length(session.epochs));
    for e = 1:length(session.epochs)
        if ~isfield(session.epochs{1,e},'startTime')
            if ~isnumeric(session.general.duration)
                sessEp{e} = Restrict(HSEuse.peaks(:), [0 str2num(session.general.duration)]);
            else
                sessEp{e} = Restrict(HSEuse.peaks(:), [0 session.general.duration]);
            end
            xtickuse{e} = 'Full session';
            warning('No stop time available');
        else
            if ~isfield(session.epochs{1,e},'stopTime')
                if ~isnumeric(session.general.duration)
                    sessEp{e} = Restrict(HSEuse.peaks(:), [session.epochs{1,e}.startTime str2num(session.general.duration)]);
                else
                    sessEp{e} = Restrict(HSEuse.peaks(:), [session.epochs{1,e}.startTime session.general.duration]);
                end
                xtickuse{e} = 'Full session';
                warning('No stop time available');
            else
                try
                    sessEp{e} = Restrict(HSEuse.peaks(:), [session.epochs{1,e}.startTime session.epochs{1,e}.stopTime]);
                    xtickuse{e} = session.epochs{1,e}.name;
                catch
                    xtickuse{e} = '';
                end
            end
        end
        bary(e) = length(sessEp{e});
    end

    figure(7);
    hold on
    bar(bary);
    xticks([1:length(session.epochs)]);
    xticklabels(xtickuse);
    title('Distribution of events across session epochs');
    ylabel('Number of events');
    labelsState = string(bary);
    text(1:length(session.epochs),bary,labelsState,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    hold off
    saveas(gcf,[plotSave 'epochSplit.png']);

end

else
    disp('no barrages during nrem');
end
end
%%%%
% function plotPSTH(PSTH, regTot, subNum, eventName,figNum)
%     [~,index3] = sort(PSTH.modulationPeakResponseTime);
%     figure(figNum);
%     subplot(4,regTot,subNum);
%     plot(PSTH.time,mean(PSTH.responsecurve')','LineWidth',2); hold on;
%     plot(PSTH.time,mean(PSTH.responsecurve')+std(PSTH.responsecurve'),'--b');hold on;
%     plot(PSTH.time,mean(PSTH.responsecurve')-std(PSTH.responsecurve'),'--b');hold on;
%     xline(0,'--k');hold on; ylabel('mod. index');
%     title(eventName)
%     figure(figNum);
%     subplot(4,regTot,subNum+regTot)
%     imagesc(PSTH.time,[1:size(PSTH.responsecurve,2)],zscore(PSTH.responsecurve(:,index3))',[-3 3]), xlabel('time'), ylabel('units');hold on;
%     xline(0,'--k');hold on;
% end
function plotPSTH(PSTH, t, index, regTot, subNum, eventName,figNum)
%     [~,index3] = sort(PSTH.modulationPeakResponseTime);
figure(figNum);
subplot(4,regTot,subNum);
plot(t,zscore(mean(PSTH,1))','LineWidth',2); hold on;
plot(t,mean(PSTH,1)+std(PSTH,1),'--b');hold on;
plot(t,mean(PSTH,1)-std(PSTH,1),'--b');hold on;
xline(0,'--k');hold on; ylabel('mod. index');
title(eventName)
figure(figNum);
subplot(4,regTot,subNum+regTot)
imagesc(t,[1:size(PSTH,1)],PSTH(index,:)); xlabel('time (s)'); ylabel('units');hold on;
% clim([-0.5 0.5]);
%     clim([min(PSTH,[],'all') max(PSTH,[],'all')]);
xline(0,'--k');hold on;
end

function [wake_theta] = getWakeStates(SleepState)
% pull all wake/theta, wake/nontheta, NREM, REM
wake_theta = []; wt = 1;
if isfield(SleepState.ints, 'WAKEstate')&&isfield(SleepState.ints, 'THETA')
    for i = 1:size(SleepState.ints.WAKEstate,1)
        tempS = []; tempS = find(SleepState.ints.THETA(:,1) >= SleepState.ints.WAKEstate(i,1));
        tempE = []; tempE = find(SleepState.ints.THETA(:,2) <= SleepState.ints.WAKEstate(i,2));
        [~,ia] = intersect(tempS, tempE);
        for j = 1:length(ia)
            wake_theta(wt,1) = SleepState.ints.THETA(tempS(ia(j)),1); wake_theta(wt,2) = SleepState.ints.THETA(tempS(ia(j)),2);
            wt = wt+1;
        end
    end
    %     if isfield(SleepState.ints, 'THETAtask')
    %         for i = 1:size(SleepState.ints.WAKEstate,1)
    %             tempS = []; tempS = find(SleepState.ints.THETAtask(:,1) >= SleepState.ints.WAKEstate(i,1));
    %             tempE = []; tempE = find(SleepState.ints.THETAtask(:,2) <= SleepState.ints.WAKEstate(i,2));
    %             [~,ia] = intersect(tempS, tempE);
    %             for j = 1:length(ia)
    %                 wake_theta(wt,1) = SleepState.ints.THETAtask(tempS(ia(j)),1); wake_theta(wt,2) = SleepState.ints.THETAtask(tempS(ia(j)),2);
    %                 wt = wt+1;
    %             end
    %         end
    %         [wake_theta(:,1), iwt] = sort(wake_theta(:,1));
    %         wake_theta(:,2) = wake_theta(iwt,2);
    %     end
end
end

function [order] = getSort(PSTH)
tempMax = max(PSTH,[],2);
[~,order] = sort(tempMax);
end

function [zscored] = manualZscore(matrix,spikes)
%matrix must be units X times dimension
zscored = matrix;
for i = 1:size(matrix,1)
    [~,spkmean,spkstd] = spkRtHist(spikes.times{i}, 'ifz',false);
    zscored(i,1:size(matrix,2)) = (matrix(i,:) - spkmean)/spkstd;
end
end

function [PSTH, PSTH_t] = manualPSTH(spikes,evtTime,win,bin)
PSTH_t = [win(1):bin:win(2)];
lenMat = length(PSTH_t);
PSTH = NaN(size(spikes.times,2), lenMat);
parfor i = 1:size(spikes.times,2)
    PSTH_cell{i} = getEvtHist(spikes.times{i}, evtTime, win, bin);
end
for i = 1:length(PSTH_cell)
    PSTH(i,1:lenMat) = PSTH_cell{i};
end
end

function evtHist = getEvtHist(spkTimes, evtTime, win, bin)
evtHist = [];
time1 = evtTime+win(1); time2 = evtTime+win(2);
tmp_spk = relativeRestrict(spkTimes,[time1 time2]);
tmp_spk = tmp_spk + win(1);
timebins = [win(1):bin:win(2)];
tmp = histc(tmp_spk,timebins)/length(tmp_spk); %divide by num spks
% tmp = tmp/length(evtTime); %divide by num events
evtHist = tmp/bin; %divide by bin size for FR output
fixIt = find(isnan(evtHist));
evtHist(fixIt)=0;
end

function spks = relativeRestrict(spkTimes, intervals)
% this should give us a list of spike times relative to each event
spks = [];
for i = 1:size(intervals,1)
    tempSpks = find(spkTimes>=intervals(i,1)&spkTimes<=intervals(i,2));
    tempSpks = spkTimes(tempSpks) - intervals(i,1);
    spks = [spks; tempSpks];
end
end