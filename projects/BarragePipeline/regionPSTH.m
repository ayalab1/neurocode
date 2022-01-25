%% Multi-PSTH
% Plot multiple PSTHs for multiple regions, barrage vs ripples

%% 
close all
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
load(strcat(basepath,'\Barrage_Files\',basename,'.HSE.mat'));
load(strcat(basename,'.ripples.events.mat'));
load(strcat(basename,'.cell_metrics.cellinfo.mat'));
savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
% plotSave = strcat('Z:\home\Lindsay\Barrage\PSTH\',animName,'.',basename,'.');
plotSave = strcat(basepath, '\Barrage_Files\', basename, '.');
%%
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
            spikes = importSpikes('brainRegion', check(i), 'cellType', "Pyramidal Cell");
            save([savePath br 'pyr.cellinfo.mat'], 'spikes');
            if ~isempty(spikes.times)
                PSTH_ripples = computePSTH(ripples,spikes,'duration',2,'plots', false);
                plotPSTH(PSTH_ripples, regTot, subNum, strcat(check(i), '/ripples'),1);
                PSTH_bar = computePSTH(HSE,spikes,'duration',2,'plots', false);
                plotPSTH(PSTH_bar, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'),1);
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
            spikes = importSpikes('brainRegion', check(i), 'cellType', ["Wide Interneuron"; "Narrow Interneuron"]);
            save([savePath br 'int.cellinfo.mat'], 'spikes');
            if ~isempty(spikes.times)
                PSTH_ripples = computePSTH(ripples,spikes,'duration',2,'plots', false);
                plotPSTH(PSTH_ripples, regTot, subNum, strcat(check(i), '/ripples'),2);
                PSTH_bar = computePSTH(HSE,spikes,'duration',2,'plots', false);
                plotPSTH(PSTH_bar, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'),2);
            else
                warning(['No interneurons in ' br]);
            end
            subNum = subNum+1; % should max at 4
        end
    end
end
sgtitle('Interneuron PSTH');
saveas(gcf,[plotSave 'PSTHint.png']);

%% Event duration
duration = HSE.timestamps(:,2)-HSE.timestamps(:,1);
figure(3);
hist(duration,50);title('barrages');xlabel('duration (s)');
saveas(gcf,[plotSave 'evtDist.png']);
%% CCG
t_ripple{1}=ripples.timestamps(:,1);
t_barrages{1}=HSE.timestamps(:,1);
t_ripple_id{1}=ones(1,length(ripples.timestamps(:,1)))';
t_barrage_id{1}=ones(1,length(HSE.timestamps(:,1)))';
binsize=0.1;duration=20;
[ccg_ripple_barrage,t_ripple_barrage] = CCG(cat(1,t_ripple{1},t_barrages{1}),cat(1,t_ripple_id{1},2*t_barrage_id{1}),'binSize',binsize,'duration',duration,'norm','rate');
figure(4);plot(t_ripple_barrage,ccg_ripple_barrage(:,2,1),'r');hold on;title('ccg barr-SWR');
saveas(gcf,[plotSave 'CCG.png']);

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