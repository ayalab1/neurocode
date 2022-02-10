%% Cumulative Metrics for Figures 

function NewMetComb(combine)
if nargin < 1
    combine = "both";
end
if ~((combine == "mouse")||(combine == "rat")||(combine == "both"))||(nargin < 1)
    warning('Must choose mouse, rat, or both. Defaulting to both');
    combine = "both";
end
savePath = ('Z:\home\Lindsay\Barrage\');
if combine == "both"
    load('Z:\home\Lindsay\Barrage\combinedPaths.mat');
elseif combine == "mouse"
    load('Z:\home\Lindsay\Barrage\mousePaths.mat');
else
    load('Z:\home\Lindsay\Barrage\ratPaths.mat');
end

NEWEVTS = 0;
%% Initialize cells and arrays
evtDur = []; %event durations
absFRpx = []; %absolute firing rate (during bar) data pyramidal cells
absFRpg = []; %absolute firing rate (during bar) group pyramidal cells
absFRix = []; %absolute firing rate (during bar) data interneurons
absFRig = []; %absolute firing rate (during bar) group interneurons
absFR2x = []; %CA2 modulation
absFR2g = []; %CA2 modulation
gainFRpb = cell(8,1);
gainFRps = cell(8,1);
gainMetp = cell(8,1);
gainMetpg = cell(8,1);
gainFRib = cell(8,1);
gainFRis = cell(8,1);
% gainMeti = cell(8,1);
burstcp = cell(8,1);

%% Iterate through saved sessions
for p = 1:size(paths_save,1) 
% for p = 1:2
%     disp(paths_save(p)); %just need some form of progress bar
%     fprintf('This message is sent at time %s\n', datestr(now,'HH:MM:SS'));
    ses = paths_save(p);
    cd(ses);
    basepath = pwd;
    basename = basenameFromBasepath(basepath);
    animName = animalFromBasepath(basepath);
    if contains(basepath, 'Z:\')
        label = "rat";
    elseif contains(basepath, 'Y:\')
        label = "mouse";
    else
        error('Check location labels');
    end
    if (label==combine)||(combine=="both")
%         load(strcat(basepath,'\Barrage_Files\Population\',basename,'.brstDt.props.mat'));
        load(strcat(basepath,'\',basename,'.cell_metrics.cellinfo.mat'));
        load(strcat(basepath,'\Barrage_Files\Cell\',basename,'.props.mat'),'regID','burstIndex','regKey');
        load(strcat(basepath,'\Barrage_Files\Population\',basename,'.brstDt.spkEventTimes.mat'));
        load(strcat(basepath,'\Barrage_Files\Population\',basename,'.brstDt.unitBar.mat'));
        load(strcat(basepath,'\',basename,'.session.mat'));
        load(strcat(basepath,'\',basename,'.spikes.cellinfo.mat')); spikesALL = spikes; clear spikes
        load(strcat(basepath,'\Barrage_Files\',basename,'.allpyr.cellinfo.mat')); spikesPYR = spikes; clear spikes
        load(strcat(basepath,'\Barrage_Files\',basename,'.allint.cellinfo.mat')); spikesINT = spikes; clear spikes
                
        intSave = strcat(basepath, '\Barrage_Files\Population\',basename,'.brstDt');
        %% Rerun mets for interneurons
        if (~exist(strcat(intSave,'.spkEventTimesINT.mat')))||(NEWEVTS)
            load(strcat(basepath,'\Barrage_Files\',basename,'.HSEfutEVT.mat'));
            [spkEventTimesINT] = getRipSpikes('basepath',pwd,'events',evtSave{end,1},'spikes',spikesINT,'padding',0,'saveMat',false);
            save(strcat(intSave,'.spkEventTimesINT.mat'),'spkEventTimesINT', '-v7.3');
            [unitBarINT] = unitSWRmetrics(spkEventTimesINT,spikesINT);
            save(strcat(intSave,'.unitBarINT.mat'),'unitBarINT', '-v7.3');
        else
            load(strcat(intSave,'.spkEventTimesINT.mat'));
            load(strcat(intSave,'.unitBarINT.mat'));
        end
        
        %% Population Level Event Duration Histogram
        evtDur = [evtDur; spkEventTimes.EventDuration];
        sesDur = (session.epochs{1,end}.stopTime-session.epochs{1,1}.startTime);
        %% PYR RUNS
        % Absolute FR during barr (per region, box)
        meanPYR = NaN(1,size(spikesPYR.UID,2));
        stdPYR = NaN(1,size(spikesPYR.UID,2));
        spkhist = cell(1,size(spikesPYR.UID,2));
        useUn = cell(1,size(spikesPYR.UID,2));
        absFRpxt = cell(1,size(spikesPYR.UID,2));
        absFRpgt = cell(1,size(spikesPYR.UID,2));
        absFR2xt = cell(1,size(spikesPYR.UID,2));
        absFR2gt = cell(1,size(spikesPYR.UID,2));
        nSpkEach = unitBar.nSpkEach; EventDuration = spkEventTimes.EventDuration;
        Narr = cell_metrics.tags.N; Parr = cell_metrics.tags.P;
        parfor u = 1:size(spikesPYR.UID,2)
            % Absolute FR during barr (per region, box)
            tSmooth = 0.02;
            binSz = 0.005;
            spkhist{u} = spkRtHist(spikesPYR.times{u}, 'tSmooth',tSmooth,'binSz',binSz,'ifz', false);
            meanPYR(u) = mean(spkhist{u}/binSz); %divide by binsz
            stdPYR(u) = std(spkhist{u}/binSz);
            useUn{u} = find(nSpkEach(u,:) ~= 0);
            if ~isempty(useUn{u})
                absFRpxt{u} = ((nSpkEach(u,useUn{u})'./EventDuration(useUn{u}))-meanPYR(u))/stdPYR(u);
                absFRpgt{u}= ones(length(useUn{u}),1)*regID(spikesPYR.UID(u));
                if regID(spikesPYR.UID(u))==3
                    absFR2xt{u} = absFRpxt{u};
                    absFR2gt{u} = sortMod(u, length(absFRpxt{u}), Narr, Parr);
                end
            end
        end
        absFRpx = [absFRpx; cat(1,absFRpxt{:})]; absFRpg = [absFRpg; cat(1,absFRpgt{:})]; 
        absFR2x = [absFR2x; cat(1,absFR2xt{:})]; absFR2g = [absFR2g; cat(1,absFR2gt{:})]; 
        clear absFRpxt absFRpgt absFR2xt absFR2gt useUn spkhist
        for u = 1:size(spikesPYR.UID,2)
            %barrage abs FR
            tempArr = []; tempArr = gainFRpb{regID(spikesPYR.UID(u))};
            tempArr = [tempArr; mean(nSpkEach(u,:)'./EventDuration(:))];
            gainFRpb{regID(spikesPYR.UID(u))} = tempArr;
            %session FR
            tempArr = []; tempArr = gainFRps{regID(spikesPYR.UID(u))};
            tempArr = [tempArr; length(spikesPYR.times{u})/sesDur];
            gainFRps{regID(spikesPYR.UID(u))} = tempArr;
            %single gain metric
            tempArr = []; tempArr = gainMetp{regID(spikesPYR.UID(u))};
            gainMet = []; 
                % (barFR - sessFR)/(barFR + sessFR)
            gainMet = ((mean(nSpkEach(u,:)'./EventDuration(:)))-(length(spikesPYR.times{u})/sesDur))/((mean(nSpkEach(u,:)'./EventDuration(:)))+(length(spikesPYR.times{u})/sesDur));
            tempArr = [tempArr; gainMet];
            gainMetp{regID(spikesPYR.UID(u))} = tempArr;
            tempArr = []; tempArr = gainMetpg{regID(spikesPYR.UID(u))};
            tempArr = [tempArr; regID(spikesPYR.UID(u))];
            gainMetpg{regID(spikesPYR.UID(u))} = tempArr;
            %burst props
            tempArr = []; tempArr = burstcp{regID(spikesPYR.UID(u))};
            tempArr = [tempArr; burstIndex(u)];
            burstcp{regID(spikesPYR.UID(u))} = tempArr;
        end 
        %% INT RUNS
%         fprintf('This message is sent at time %s\n', datestr(now,'HH:MM:SS'));
        % Absolute FR during barr (per region, box)
        meanINT = NaN(1,size(spikesINT.UID,2));
        stdINT = NaN(1,size(spikesINT.UID,2));
        spkhist = cell(1,size(spikesINT.UID,2));
        useUn = cell(1,size(spikesINT.UID,2));
        absFRixt = cell(1,size(spikesINT.UID,2));
        absFRigt = cell(1,size(spikesINT.UID,2));
        nSpkEach = unitBarINT.nSpkEach; EventDuration = spkEventTimesINT.EventDuration;
        parfor u = 1:size(spikesINT.UID,2)
            % Absolute FR during barr (per region, box)
            tSmooth = 0.02;
            binSz = 0.005;
            spkhist{u} = spkRtHist(spikesINT.times{u}, 'tSmooth',tSmooth,'binSz',binSz,'ifz', false);
            meanINT(u) = mean(spkhist{u}/binSz); %divide by binsz
            stdINT(u) = std(spkhist{u}/binSz);
            % session FR
            useUn{u} = find(nSpkEach(u,:) ~= 0);
            if ~isempty(useUn{u})
                absFRixt{u} = ((nSpkEach(u,useUn{u})'./EventDuration(useUn{u}))-meanINT(u))/stdINT(u);
                absFRigt{u}= ones(length(useUn{u}),1)*regID(spikesINT.UID(u));
            end
        end
        absFRix = [absFRix; cat(1,absFRixt{:})]; absFRig = [absFRig; cat(1,absFRigt{:})]; clear absFRixt absFRigt useUn spkhist
        for u = 1:size(spikesINT.UID,2)
            %barrage abs FR
            tempArr = []; tempArr = gainFRib{regID(spikesINT.UID(u))};
            tempArr = [tempArr; mean(nSpkEach(u,:)'./EventDuration(:))];
            gainFRib{regID(spikesINT.UID(u))} = tempArr;
            %session FR
            tempArr = []; tempArr = gainFRis{regID(spikesINT.UID(u))};
            tempArr = [tempArr; length(spikesINT.times{u})/sesDur];
            gainFRis{regID(spikesINT.UID(u))} = tempArr;
        end
    end 
end

regKeyUse = [];
for i = 1:size(regKey,2)
    regKeyUse(i) = str2double(regKey(2,i));
end

%% %%% Plotting %%% %%

%% Event duration histogram
figure('Position', get(0, 'Screensize'));
binsEvtDur = [0:0.01:2];
evtDurHist = histc(evtDur, binsEvtDur);
% evtDurHistSmooth = smooth(evtDurHist, 'sgolay', 1);
% plot(binsEvtDur, evtDurHist./sum(evtDurHist));
plot(binsEvtDur, evtDurHist);
xlim([0 2]);
% plot(binsEvtDur, evtDurHist);
title('Barrage Event Duration');
ylabel('Barrage Count');
xlabel('Event duration (s)');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.evtDurHist.png']);

%% PYR Absolute FR during barrages (boxplot)
figure('Position', get(0, 'Screensize'));
boxplot(absFRpx, absFRpg);
[~,useLabels] = intersect(regKeyUse,absFRpg); 
xticklabels(regKey(1,useLabels));
zline = refline(0,0);
zline.LineStyle = '--'; zline.Color = 'k';
title('Pyramidal Cell absolute firing rate during barrages, per region');
ylabel('Absolute Firing Rate (# spikes/second)');
xlabel('Region');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRabsFR.png']);

%% INT Absolute FR during barrages (boxplot)
figure('Position', get(0, 'Screensize'));
boxplot(absFRix, absFRig);
[~,useLabels] = intersect(regKeyUse,absFRig); 
xticklabels(regKey(1,useLabels));
zline = refline(0,0);
zline.LineStyle = '--'; zline.Color = 'k';
title('Interneuron absolute firing rate during barrages, per region');
ylabel('Absolute Firing Rate (# spikes/second)');
xlabel('Region');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.INTabsFR.png']);

%% CA2Mod Absolute FR during barrages (boxplot)
figure('Position', get(0, 'Screensize'));
boxplot(absFR2x, absFR2g);
modLabels = ["Unknown"; "P"; "N"];
xticklabels(modLabels(unique(absFR2g)));
zline = refline(0,0);
zline.LineStyle = '--'; zline.Color = 'k';
title('Pyramidal CA2 absolute firing rate during barrages, per modulation type');
ylabel('Absolute Firing Rate (# spikes/second)');
xlabel('CA2 Modulation Type');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.CA2absFR.png']);

%% PYR FR Gain
figure('Position', get(0, 'Screensize'));
legendUse = [];
hold on
for i = 1:length(regKeyUse)
    if ~isempty(gainFRpb{regKeyUse(i)})
        scatter(gainFRpb{regKeyUse(i)}, gainFRps{regKeyUse(i)});
%         plot(10E-3:10E-3:50, log10(10E-3:10E-3:50), 'Color','k');
        set(gca,'xscale','log', 'yscale','log');
%         set(gca,'xscale','log');
        ylim([0 20]);
        xlim([0 100]);
        legendUse = [legendUse regKey(1,i)];
    end
end
% unity=refline(1,0);
% unity.Color = 'k';
title('Firing Rate Gain PYR');
xlabel('Absolute FR during Barrage');
ylabel('Average FR across session');
legend(legendUse);
plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
hold off
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainFR.png']);

%subplots
figure('Position', get(0, 'Screensize'));
for i = 1:length(regKeyUse)
    subplot(4,2,i);
    hold on
    title([regKey(1,i)]);
    if ~isempty(gainFRpb{regKeyUse(i)})
        scatter(gainFRpb{regKeyUse(i)}, gainFRps{regKeyUse(i)});
        set(gca,'xscale','log', 'yscale','log');
%         unity=refline(1,0);
%         unity.Color = 'k';
        xlabel('Absolute FR during Barrage');
        ylabel('Average FR across session');
    end
    plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
    ylim([10E-4 20]);
    xlim([10E-4 50]);
    hold off
end
sgtitle('Firing Rate Gain PYR');

saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainFRsub.png']);

% single metric per region box
gainMetp = cat(1,gainMetp{:}); gainMetpg = cat(1,gainMetpg{:});
figure('Position', get(0, 'Screensize'));
boxplot(gainMetp, gainMetpg);
yline(0, '--');
[~,useLabels] = intersect(regKeyUse,gainMetpg); 
xticklabels(regKey(1,useLabels));
title('Gain Metric per region');
ylabel('(BarrFR-SessFR)/(BarrFR+SessFR)');
xlabel('Region');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainMetBox.png']);

%% Interneuron 
figure('Position', get(0, 'Screensize'));
legendUse = [];
hold on
for i = 1:length(regKeyUse)
    if ~isempty(gainFRib{regKeyUse(i)})
        scatter(gainFRib{regKeyUse(i)}, gainFRis{regKeyUse(i)});
        set(gca,'xscale','log', 'yscale','log');
        legendUse = [legendUse regKey(1,i)];
    end
end
% unity=refline(1,0);
% unity.Color = 'k';
title('Firing Rate Gain INT');
xlabel('Absolute FR during Barrage');
ylabel('Average FR across session');
legend(legendUse);
plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
% ylim([0 15]);
% xlim([0 100]);
hold off
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.INTgainFR.png']);

%subplots
figure('Position', get(0, 'Screensize'));
for i = 1:length(regKeyUse)
    subplot(4,2,i);
    hold on
    title([regKey(1,i)]);
    if ~isempty(gainFRib{regKeyUse(i)})
        scatter(gainFRib{regKeyUse(i)}, gainFRis{regKeyUse(i)});
        set(gca,'xscale','log', 'yscale','log');
%         unity=refline(1,0);
%         unity.Color = 'k';
        xlabel('Absolute FR during Barrage');
        ylabel('Average FR across session');
    end
    plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
%     ylim([0 15]);
%     xlim([0 100]);
    hold off
end

sgtitle('Firing Rate Gain INT');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.INTgainFRsub.png']);

%% Burst Properties PYR ONLY
figure('Position', get(0, 'Screensize'));
for i = 1:length(regKeyUse)
    subplot(4,2,i);
    hold on
    title([regKey(1,i)]);
    if ~isempty(burstcp{regKeyUse(i)})
        scatter(gainFRpb{regKeyUse(i)}, burstcp{regKeyUse(i)});
        set(gca,'xscale','log', 'yscale','log');
        xlabel('Absolute FR during Barrage');
        ylabel('Burst Index');
        xlim([10E-4 10]); 
        ylim([10E-4 1]);
    end
    hold off
end
sgtitle('Absolute FR during Barrage vs Burst Index of Cells');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRcellBurstsub.png']);

close all

end


function group = sortMod(u, lenGroup, Narr, Parr)
    if ~isempty(find(Narr==u,1))
        group = ones(lenGroup,1)*3;
    elseif ~isempty(find(Parr==u,1))
        group = ones(lenGroup,1)*2;
    else
        group = ones(lenGroup,1)*1;
    end
end