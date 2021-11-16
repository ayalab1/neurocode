%% Multi-PSTH
% Plot multiple PSTHs for multiple regions, barrage vs ripples

%% 
clear
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
load(strcat(basepath,'\Barrage_Files\',basename,'.HSE.mat'));
load(strcat(basename,'.ripples.events.mat'));
load(strcat(basename,'.cell_metrics.cellinfo.mat'));
savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
%%
regions = unique(cell_metrics.brainRegion);
check = ["CA1" "CA2" "CA3" "CTX" "DG"];
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
            PSTH_ripples = computePSTH(ripples,spikes,'duration',2,'plots', false);
            plotPSTH(PSTH_ripples, regTot, subNum, strcat(check(i), '/ripples'));
            PSTH_bar = computePSTH(HSE,spikes,'duration',2,'plots', false);
            plotPSTH(PSTH_bar, regTot, subNum+(regTot*2), strcat(check(i), '/barrages'));
            subNum = subNum+1; % should max at 4
        end
    end
end
%%%%
function plotPSTH(PSTH, regTot, subNum, eventName)
    [~,index3] = sort(PSTH.modulationPeakResponseTime);  
    figure(1);
    subplot(4,regTot,subNum);
    plot(PSTH.time,mean(PSTH.responsecurve')','LineWidth',2); hold on; 
    plot(PSTH.time,mean(PSTH.responsecurve')+std(PSTH.responsecurve'),'--b');hold on; 
    plot(PSTH.time,mean(PSTH.responsecurve')-std(PSTH.responsecurve'),'--b');hold on;
    xline(0,'--k');hold on; ylabel('mod. index');
    title(eventName)
    figure(1);
    subplot(4,regTot,subNum+regTot)
    imagesc(PSTH.time,[1:size(PSTH.responsecurve,2)],zscore(PSTH.responsecurve(:,index3))',[-3 3]), xlabel('time'), ylabel('units');hold on;
    xline(0,'--k');hold on;    
end