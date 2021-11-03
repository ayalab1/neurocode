%% Cumulative Metrics

% We're essentially just going to pull the box plot data and go from there

% check if input is a mouse or a rat file (contains 'Y:\' or 'Z:\' in
% basepath for session)

% Let's do this separately from BarAnalysis and use the session list to
% remake as we go. Data is already there, just concatenating and plotting
function CumMetRun(combine)
if ~((combine == "mouse")||(combine == "rat")||(combine == "both"))
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
avgSpk_boxx = [];
avgSpk_boxg = [];
avgFR_boxx = [];
avgFR_boxg = [];

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
        load(strcat(basepath,'\Barrage_Files\',basename,'.allpyr.cellinfo.mat'));
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
        burstSz_boxx = [burstSz_boxx; cumMet.boxxBurstSz];
        burstSz_boxg = [burstSz_boxg; cumMet.boxgBurstSz];
        
        %% Cell level histogram of burst spike lengths
        for r = 1:size(burstLen,1)
            for c = 1:size(burstLen,2)
                if (size(cumMet.burstCnt,1)>=r)&&(size(cumMet.burstCnt,2)>=c)
                    burstLen{r,c} = [burstLen{r,c} cumMet.burstCnt{r,c}];
                end
            end
        end
        
        %% Population Level Percent of Units in Events
        perc_boxx = [perc_boxx cumMet.boxxPerc];
        perc_boxg = [perc_boxg cumMet.boxgPerc];

        %% Population Level Avg Spike per Unit
        avgSpk_boxx = [avgSpk_boxx; cumMet.boxxAvgSpk];
        avgSpk_boxg = [avgSpk_boxg; cumMet.boxgAvgSpk];

        %% Population Level Avg FR per Unit
        avgFR_boxx = [avgFR_boxx; cumMet.boxxAvgFR];
        avgFR_boxg = [avgFR_boxg; cumMet.boxgAvgFR];
    end 
end

check = ["CA1" "CA2" "CA3" "CTX" "DG"];
useName = unique(ISI_boxg);
nameInd = ismember(1:length(check),useName);
presentIND = find(nameInd==1);
presentName = convertStringsToChars(check(nameInd));

%% Plotting 

% Cell level ISI dist
% [ISI,ISIc,t] = ISIGrams(spikes.times, spikes.UID, 1/1000, 1000);
% cellProp.ISI = ISI;
% cellProp.ISIhistT = t;
% cellProp.ISIhist = ISIc;
% ISIavg = NaN(length(ISI),1);
% for i = 1:length(ISI)
%     ISIavg(i) = mean(ISI{i});
% end
% cellProp.ISIavg = ISIavg;


figure('Position', get(0, 'Screensize'));
plot(ISI_x,ISI_y);
xlabel('Log of ISI (ms)');
ylabel('Count');
title('Distribution of ISIs in Log scale');
set(gca, 'XScale', 'log')
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.ISIdistCum.png']);

% Cell level type ISI log distributions
figure('Position', get(0, 'Screensize'));
title('Log(ISI) of per region types');
c = 2;
i = 1;
for r = 1:(size(ISI_c,1))
    subplot(size(ISI_c,1),2,i); plot(ISI_x, ISI_c{r,c});hold on;title(strcat(cumMet.regKey(1,r),' P'));
    i=i+1;
    subplot(size(ISI_c,1),2,i); plot(ISI_x, ISI_c{r,c+1});hold on; title(strcat(cumMet.regKey(1,r),' N'));
    i=i+1;
end
hold off;
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.ISItypeDistCum.png']);

% Cell level ISI
figure('Position', get(0, 'Screensize'));
boxplot(ISI_boxx, ISI_boxg);
xticklabels(presentName);
title('ISI per region, outliers cut off');
ylabel('ISI (s)');
xlabel('Region');
ylim([-1 7]);
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.ISIboxCum.png']);

% Cell Level Burst Index
figure('Position', get(0, 'Screensize'));
boxplot(burst_boxx, burst_boxg);
xticklabels(presentName);
title('Burst Index per region');
ylabel('Burst Index');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.BurstBoxCum.png']);

% Cell level Burst Size
figure('Position', get(0, 'Screensize'));
boxplot(burstSz_boxx, burstSz_boxg);
xticklabels(presentName);
title('Burst Size per region');
ylabel('Burst Size (# spikes)');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.BurstSzBoxCum.png']);

% Cell level burst length histogram
figure('Position', get(0, 'Screensize'));
title('Histograms of burst length per region types');
c = 2; %we don't care about unknown modulations
i = 1;
for r = 1:(size(burstLen,1))
    subplot(size(burstLen,1),2,i); hist(burstLen{r,c},[2:2:10])/sum(hist(burstLen{r,c},[2:2:10]));hold on; title(strcat(cumMet.regKey(1,r),' P'));
    i=i+1;
    subplot(size(burstLen,1),2,i); hist(burstLen{r,c+1},[2:2:10]);hold on; title(strcat(cumMet.regKey(1,r),' N'));
    i=i+1;
end
hold off;
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.BurstHistCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(perc_boxx, perc_boxg);
xticklabels(presentName);
title('Percent of units making up event per region');
ylabel('Percent of units');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.PercUntsCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(avgSpk_boxx, avgSpk_boxg);
xticklabels(presentName);
title('Average number of spikes per event per region');
ylabel('Number of spikes');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.AvgSpkCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(avgFR_boxx, avgFR_boxg);
xticklabels(presentName);
title('Average Firing Rate per event per region');
ylabel('Firing Rate (Hz)');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(combine) '.AvgFRCum.png']);

close all

end