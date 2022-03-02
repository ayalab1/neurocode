%% Multi-PSTH
% Plot multiple PSTHs for multiple regions, barrage vs ripples

%% 
close all
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
load(strcat(basepath,'\Barrage_Files\',basename,'.HSE.mat'));

load(strcat(basename,'.cell_metrics.cellinfo.mat'));
load(strcat(basepath,'\',basename,'.SleepState.states.mat'));
savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
% plotSave = strcat('Z:\home\Lindsay\Barrage\PSTH\',animName,'.',basename,'.');
plotSave = strcat(basepath, '\Barrage_Files\', basename, '.');

if ~exist(strcat(basepath,'\Barrage_Files\PSTHmet\'))
        mkdir(strcat(basepath,'\Barrage_Files\PSTHmet\'));
end
if ~exist(strcat(basepath,'\',basename,'.ripples.events.mat'))
    if exist(strcat(basepath,'\old files\SWR.mat'))
        load(strcat(basepath,'\old files\SWR.mat'));
        ripples = SWR;
        if isfield(ripples,'peaktimes')
            ripples.peaks = ripples.peaktimes;
        end
        plotRips = 1;
    elseif exist(strcat(basepath,'\oldfiles\SWR.mat'))
        load(strcat(basepath,'\oldfiles\SWR.mat'));
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
    plotRips = 1;
end
%%

if isfield(HSE, 'keep')
    HSE.timestamps = HSE.timestamps(HSE.keep,:);
    HSE.peaks = HSE.peaks(HSE.keep);
    HSE.amplitudes = HSE.amplitudes(HSE.keep);
    HSE.duration = HSE.duration(HSE.keep);
end
HSEnrem = eventIntervals(HSE,SleepState.ints.NREMstate,1);

regions = unique(cell_metrics.brainRegion);
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

subNum = 1;
for i = 1:length(check)
    for j=1:length(regions)
        regCheck = regions{j};
        if contains(regCheck, check(i))
            br = convertStringsToChars(check(i));
            if ~exist([savePath br 'pyr.cellinfo.mat'])
                spikes = importSpikes('brainRegion', check(i), 'cellType', "Pyramidal Cell");
                save([savePath br 'pyr.cellinfo.mat'], 'spikes');
            else 
                load([savePath br 'pyr.cellinfo.mat']);
            end
            if ~isempty(spikes.times)
                if plotRips
                    PSTH_ripples = computePSTH(ripples,spikes,'duration',2,'plots', false,'alignment','peaks');
                    plotPSTH(PSTH_ripples, regTot, subNum, strcat(check(i), '/ripples'),1);
                    PSTHmets.PSTH_ripples = PSTH_ripples;
                end
                    PSTH_bar = computePSTH(HSEnrem,spikes,'duration',2,'plots', false,'alignment','peaks');
                    plotPSTH(PSTH_bar, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'),1);
                    PSTHmets.PSTH_bar = PSTH_bar;
                    save([basepath '\Barrage_Files\PSTHmet\' basename '.' br 'PSTHmets.mat'], 'PSTHmets');
                else
                    warning(['No pyramidal cells in ' br]);
                end
            
            subNum = subNum+1; % should max at 4
        end
    end
end
sgtitle('Pyramidal PSTH');
saveas(gcf,[plotSave 'PSTH.png']);

subNum = 1;
for i = 1:length(check)
    for j=1:length(regions)
        regCheck = regions{j};
        if contains(regCheck, check(i))
            br = convertStringsToChars(check(i));
            spikes = [];
            if ~exist([savePath br 'int.cellinfo.mat'])
                spikes = importSpikes('brainRegion', check(i), 'cellType', ["Wide Interneuron"; "Narrow Interneuron"]);
                save([savePath br 'int.cellinfo.mat'], 'spikes');
            else
                load([savePath br 'int.cellinfo.mat']);
            end
            if ~isempty(spikes.times)
                if plotRips
                    PSTH_ripples = computePSTH(ripples,spikes,'duration',2,'plots', false,'alignment','peaks');
                    plotPSTH(PSTH_ripples, regTot, subNum, strcat(check(i), '/ripples'),2);
                    PSTHmetsINT.PSTH_ripples = PSTH_ripples;
                end
                PSTH_bar = computePSTH(HSEnrem,spikes,'duration',2,'plots', false,'alignment','peaks');
                plotPSTH(PSTH_bar, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'),2);
                PSTHmetsINT.PSTH_bar = PSTH_bar;
                save([basepath '\Barrage_Files\PSTHmet\' basename '.' br 'PSTHmetsINT.mat'], 'PSTHmetsINT');
            else
                warning(['No interneurons in ' br]);
            end
            subNum = subNum+1; % should max at 4
        end
    end
end
if ~isempty(spikes.UID)
    sgtitle('Interneuron PSTH');
    saveas(gcf,[plotSave 'PSTHint.png']);
end

%% Event duration
duration = HSEnrem.timestamps(:,2)-HSEnrem.timestamps(:,1);
figure(3);
hist(duration,50);title('barrages');xlabel('duration (s)');
saveas(gcf,[plotSave 'evtDist.png']);
%% CCG
if plotRips
    t_ripple{1}=ripples.peaks(:);
    t_barrages{1}=HSEnrem.peaks(:);
    t_ripple_id{1}=ones(1,length(ripples.peaks(:)))';
    t_barrage_id{1}=ones(1,length(HSEnrem.peaks(:)))';
    binsize=0.1;duration=6;
    [ccg_ripple_barrage,t_ripple_barrage] = CCG(cat(1,t_ripple{1},t_barrages{1}),cat(1,t_ripple_id{1},2*t_barrage_id{1}),'binSize',binsize,'duration',duration,'norm','rate');
    figure(4);plot(t_ripple_barrage,ccg_ripple_barrage(:,2,1),'r');hold on;title('ccg barr-SWR');
    saveas(gcf,[plotSave 'CCG.png']);
    CCG_dat.time = t_ripple_barrage; CCG_dat.y = ccg_ripple_barrage(:,2,1);
    save([savePath 'CCG_dat.mat'],'CCG_dat');
end

%% Fraction per state
[wake_theta] = getWakeStates(SleepState);
num_theta = 0; num_nontheta = 0; num_REM = 0; num_NREM = 0;
if ~isempty(wake_theta)
    temp_theta = Restrict(HSE.peaks(:),wake_theta);
else
    temp_theta = [];
end
num_theta = length(temp_theta);
temp_wake = Restrict(HSE.peaks(:),SleepState.ints.WAKEstate);
num_nontheta = length(temp_wake) - num_theta; %if it's wake but not theta, must be nontheta wake
num_REM = length(Restrict(HSE.peaks(:),SleepState.ints.REMstate));
num_NREM = length(Restrict(HSE.peaks(:),SleepState.ints.NREMstate));
figure(5); bar([num_theta, num_nontheta, num_REM, num_NREM]);
xticklabels({'Theta', 'Nontheta', 'REM', 'NREM'});
title('Distribution of events across states');
ylabel('Number of events');
labelsState = string([num_theta, num_nontheta, num_REM, num_NREM]);
text([1 2 3 4],[num_theta, num_nontheta, num_REM, num_NREM],labelsState,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
saveas(gcf,[plotSave 'stateSplit.png']);    

%% Fraction per task state
% maybe we could change this to be a stacked bar plot per task portion (so
% pre on x, then stack wake/nontheta/rem/nrem
load([basename '.session.mat']);
sessEp = cell(1, length(session.epochs));
xtickuse = cell(1, length(session.epochs));
bary = nan(1, length(session.epochs));
for e = 1:length(session.epochs)
    if ~isfield(session.epochs{1,e},'startTime') 
        if ~isnumeric(session.general.duration)
            sessEp{e} = Restrict(HSEnrem.peaks(:), [0 str2num(session.general.duration)]);
        else
            sessEp{e} = Restrict(HSEnrem.peaks(:), [0 session.general.duration]);
        end
        xtickuse{e} = 'Full session';
        warning('No stop time available');
    else
        if ~isfield(session.epochs{1,e},'stopTime')
            if ~isnumeric(session.general.duration)
                sessEp{e} = Restrict(HSEnrem.peaks(:), [session.epochs{1,e}.startTime str2num(session.general.duration)]);
            else
                sessEp{e} = Restrict(HSEnrem.peaks(:), [session.epochs{1,e}.startTime session.general.duration]);
            end
            xtickuse{e} = 'Full session';
            warning('No stop time available');
        else
            sessEp{e} = Restrict(HSEnrem.peaks(:), [session.epochs{1,e}.startTime session.epochs{1,e}.stopTime]);  
            xtickuse{e} = session.epochs{1,e}.name;
        end
    end
    bary(e) = length(sessEp{e});
end

figure(6);
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

%%%%
function plotPSTH(PSTH, regTot, subNum, eventName,figNum)
    [~,index3] = sort(PSTH.modulationPeakResponseTime);  
    figure(figNum);
    subplot(4,regTot,subNum);
    plot(PSTH.time,mean(PSTH.responsecurve')','LineWidth',2); hold on; 
    plot(PSTH.time,mean(PSTH.responsecurve')+std(PSTH.responsecurve'),'--b');hold on; 
    plot(PSTH.time,mean(PSTH.responsecurve')-std(PSTH.responsecurve'),'--b');hold on;
    xline(0,'--k');hold on; ylabel('mod. index');
    title(eventName)
    figure(figNum);
    subplot(4,regTot,subNum+regTot)
    imagesc(PSTH.time,[1:size(PSTH.responsecurve,2)],zscore(PSTH.responsecurve(:,index3))',[-3 3]), xlabel('time'), ylabel('units');hold on;
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