%% Cumulative Metrics

% We're essentially just going to pull the box plot data and go from there

% check if input is a mouse or a rat file (contains 'Y:\' or 'Z:\' in
% basepath for session)

% Let's do this separately from BarAnalysis and use the session list to
% remake as we go. Data is already there, just concatenating and plotting
function CumMetRun(combine)
if nargin < 1
    combine = "both";
end
if ~((combine == "mouse")||(combine == "rat")||(combine == "both"))||(nargin < 1)
    warning('Must choose mouse, rat, or both. Defaulting to both');
    combine = "both";
end
savePath = ('Z:\home\Lindsay\Barrage\');
load('Z:\home\Lindsay\Barrage\combinedPaths.mat');

ISI_x = [];
ISI_y = [];
ISI_c = cell(6,3);
for i = 1:size(ISI_c,1)
    for j = 1:size(ISI_c,2)
        ISI_c{i,j} = zeros(1000,1);
    end
end
ISI_boxx = [];
ISI_boxg = [];
burst_boxx = [];
burst_boxg = [];
burstSz_boxx = [];
burstSz_boxg = [];
burstLen = cell(6,3);
perc_boxx = [];
perc_boxg = [];
spkSC = cell(6,2);
FRsc = cell(6,2);

for p = 1:size(paths_save,1)
% for p = 1:2
    ses = paths_save(p); %can change this to an iteratable variable for mass runs
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
        load(strcat('Z:\home\Lindsay\Barrage\CumMet\',animName,'.',basename,'.cumMet.mat'));
        load(strcat(basepath,'\Barrage_Files\',basename,'.brstDt.cellinfo.mat'));
        %% Cell level ISI distribution
        ISI_x = [cumMet.ISIx]; %this should be the same for everything?
        ISI_y = sum([ISI_y,cumMet.ISIy],2);
        
        %% Cell level ISI distribution per type
        for r = 1:size(ISI_c,1)
            for c = 1:size(ISI_c,2)
                reg_id = find(cumMet.regID==r);
                mod_id = find(cumMet.modID==c);
                [~,~,useIt] = intersect(reg_id,mod_id);
                for m = 1:length(useIt)
                    useInd = find(spikes.UID==mod_id(useIt(m)),1);
                    if ~isempty(find(spikes.UID==mod_id(useIt(m)),1))
                        ISI_c{r,c} = ISI_c{r,c} + cumMet.ISIc(:,useInd);
                    end
                end
                
            end
        end
        
        %% Cell level ISI box plot
        ISI_boxx = [ISI_boxx; cumMet.boxxISI];
        ISI_boxg = [ISI_boxg; cumMet.boxgISI];

        %% Cell level burst index
        burst_boxx = [burst_boxx; cumMet.boxxBurst];
        burst_boxg = [burst_boxg; cumMet.boxgBurst];
        
        %% Cell level burst spike count
        burstSz_boxx = [burstSz_boxx; cumMet.boxxAvgBurstSz];
        burstSz_boxg = [burstSz_boxg; cumMet.boxgAvgBurstSz];
        
        %% Cell level histogram of burst spike lengths
        for r = 1:size(burstLen,1)
            for c = 1:size(burstLen,2)
                if (size(cumMet.burstCnt,1)>=r)&&(size(cumMet.burstCnt,2)>=c)
                    if ~isempty(cumMet.burstCnt{r,c})
                        tempCat = [];
                        tempCat = cat(2,cumMet.burstCnt{r,c}{:});
                        burstLen{r,c} = [burstLen{r,c} tempCat];
                    end
                end
            end
        end
        
        %% Population Level Percent of Units in Events
        perc_boxx = [perc_boxx; cumMet.boxxPerc];
        perc_boxg = [perc_boxg; cumMet.boxgPerc];

        %% Population Level Avg Spike per Unit
        for r = 1:size(spkSC,1)
            for m = 1:size(spkSC,2)
                spkSC{r,m} = [spkSC{r,m}; cumMet.cumSpk{r,m+1}];
            end
        end

        %% Population Level Avg FR per Unit
        for r = 1:size(FRsc,1)
            for m = 1:size(FRsc,2)
                FRsc{r,m} = [FRsc{r,m}; cumMet.cumFR{r,m+1}];
            end
        end
    end 
end

%% Plotting 

figure('Position', get(0, 'Screensize'));
plot(ISI_x,ISI_y);
xlabel('Log of ISI (ms)');
ylabel('Count');
title('Distribution of ISIs in Log scale');
set(gca, 'XScale', 'log')
xlim([min(ISI_x) 10E3]);
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.ISIdistCum.png']);

% Cell level type ISI log distributions
figure('Position', get(0, 'Screensize'));
title('Log(ISI) of per region types');
hold on;
c = 2;
i = 1;
for r = 1:(size(ISI_c,1))
    subplot(size(ISI_c,1),2,i); plot(ISI_x, ISI_c{r,c});hold on;title(strcat(cumMet.regKey(1,r),' P'));
    xlim([0 500]);
    set(gca, 'XScale', 'log')
    i=i+1;
    subplot(size(ISI_c,1),2,i); plot(ISI_x, ISI_c{r,c+1});hold on; title(strcat(cumMet.regKey(1,r),' N'));
    xlim([0 500]);
    set(gca, 'XScale', 'log')
    i=i+1;
end
hold off;
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.ISItypeDistCum.png']);

% Cell level ISI
figure('Position', get(0, 'Screensize'));
boxplot(ISI_boxx, ISI_boxg);
useIt = [];
setBox = unique(ISI_boxg);
for i = 1:length(setBox)
    useIt = [useIt cumMet.regKey(1,setBox(i))];
end
xticklabels(useIt);
title('ISI per region, outliers cut off');
ylabel('ISI (s)');
xlabel('Region');
ylim([-1 7]);
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.ISIboxCum.png']);

% Cell Level Burst Index
figure('Position', get(0, 'Screensize'));
boxplot(burst_boxx, burst_boxg);
useIt = [];
setBox = unique(burst_boxg);
for i = 1:length(setBox)
    useIt = [useIt cumMet.regKey(1,setBox(i))];
end
xticklabels(useIt);
title('Burst Index per region');
ylabel('Burst Index');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.BurstBoxCum.png']);

% Cell level Burst Size
figure('Position', get(0, 'Screensize'));
boxplot(burstSz_boxx, burstSz_boxg);
useIt = [];
setBox = unique(burstSz_boxg);
for i = 1:length(setBox)
    useIt = [useIt cumMet.regKey(1,setBox(i))];
end
xticklabels(useIt);
title('Burst Size per region');
ylabel('Burst Size (# spikes)');
xlabel('Region');
ylim([0 10]);
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.BurstSzBoxCum.png']);

% Cell level burst length histogram
figure('Position', get(0, 'Screensize'));
title('Histograms of burst length per region types');
c = 2; %we don't care about unknown modulations
i = 1;
for r = 1:(size(burstLen,1))
    if ~isempty(burstLen{r,c})
        tempCat = [];
        tempCat = cat(2, burstLen{r,c}(:));
        subplot(size(burstLen,1),2,i); histogram(tempCat,[2:2:10]);hold on; title(strcat(cumMet.regKey(1,r),' P'));
    end
    i=i+1;
    if ~isempty(burstLen{r,c+1})
        tempCat = [];
        tempCat = cat(2, burstLen{r,c+1}(:));
        subplot(size(burstLen,1),2,i); histogram(tempCat,[2:2:10]);hold on; title(strcat(cumMet.regKey(1,r),' N'));
    end
    i=i+1;
end
hold off;
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.BurstHistCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(perc_boxx, perc_boxg);
useIt = [];
setBox = unique(perc_boxg);
for i = 1:length(setBox)
    useIt = [useIt cumMet.regKey(1,setBox(i))];
end
xticklabels(useIt);
title('Percent of units making up event per region');
ylabel('Percent of units');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.PercUntsCum.png']);

% We're gonna make this a histogram rather than sc (maybe later we can
% sort by average intensity to get a nice trend or something?
figure('Position', get(0, 'Screensize'));
i = 1;
for r = 1:size(spkSC,1)
    if ~isempty(spkSC{r,1})
        tempSpkSC = sum(spkSC{r,1},1)/sum(spkSC{r,1});
        subplot(size(spkSC,1),2,i); bar(tempSpkSC);hold on; title(strcat(cumMet.regKey(1,r),' ',cumMet.modKey(1,2)));
        xticks([1 2 3 4 5]);
        xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)});
    end
    i = i+1;
    if r == size(spkSC,1)
        ylabel('Normalized Count');
        xlabel('# spikes');
    end
    if ~isempty(spkSC{r,2})
        tempSpkSC = sum(spkSC{r,2},1)/sum(spkSC{r,2});
        subplot(size(spkSC,1),2,i);bar(tempSpkSC);hold on; title(strcat(cumMet.regKey(1,r),' ',cumMet.modKey(1,3)));
        xticks([1 2 3 4 5]);
        xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)});
    end
    i = i+1;
end
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.SpkHistCum.png']);

%
figure('Position', get(0, 'Screensize'));
i = 1;
for r = 1:size(FRsc,1)
    if ~isempty(FRsc{r,1})
        tempFRsc = sum(FRsc{r,1},1)/sum(FRsc{r,1});
        subplot(size(FRsc,1),2,i);bar(tempFRsc);hold on; title(strcat(cumMet.regKey(1,r),' ',cumMet.modKey(1,2)));
        xticks([1:10]);
        xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)...
            num2str(12) num2str(14) num2str(16) num2str(18) num2str(20)});
    end
    i = i+1;
    if r == size(FRsc,1)
        ylabel('Normalized Count');
        xlabel('FR');
    end
    if ~isempty(FRsc{r,2})
        tempFRsc = sum(FRsc{r,2},1)/sum(FRsc{r,2});
        subplot(size(FRsc,1),2,i);bar(tempFRsc);hold on; title(strcat(cumMet.regKey(1,r),' ',cumMet.modKey(1,3)));
        xticks([1:10]);
        xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)...
            num2str(12) num2str(14) num2str(16) num2str(18) num2str(20)});
    end
    i = i+1;
end
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.FRHistCum.png']);

close all

end