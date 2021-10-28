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
ISI_boxx = [];
ISI_boxg = [];
burst_boxx = [];
burst_boxg = [];
perc_boxx = [];
perc_boxg = [];
avgSpk_boxx = [];
avgSpk_boxg = [];
avgFR_boxx = [];
avgFR_boxg = [];

% for p = 1:length(paths)
for p = 16
    ses = paths(p); %can change this to an iteratable variable for mass runs
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

        %% Cell level ISI box plot
        ISI_boxx = [ISI_boxx cumMet.boxxISI];
        ISI_boxg = [ISI_boxg cumMet.boxgISI];

        %% Cell level burst index
        burst_boxx = [burst_boxx cumMet.boxxBurst];
        burst_boxg = [burst_boxg cumMet.boxgBurst];
        
        %% Population Level Percent of Units in Events
        perc_boxx = [perc_boxx cumMet.boxxPerc];
        perc_boxg = [perc_boxg cumMet.boxgPerc];

        %% Population Level Avg Spike per Unit
        avgSpk_boxx = [avgSpk_boxx cumMet.boxxAvgSpk];
        avgSpk_boxg = [avgSpk_boxg cumMet.boxgAvgSpk];

        %% Population Level Avg FR per Unit
        avgFR_boxx = [avgFR_boxx cumMet.boxxAvgFR];
        avgFR_boxg = [avgFR_boxg cumMet.boxgAvgFR];
    end 
end

check = ["CA1" "CA2" "CA3" "CTX" "DG"];
useName = unique(ISI_boxg);
nameInd = ismember(1:length(check),useName);
presentIND = find(nameInd==1);
presentName = convertStringsToChars(check(nameInd));

%% Plotting 
% Cell level ISI
figure('Position', get(0, 'Screensize'));
boxplot(ISI_boxx, ISI_boxg);
xticklabels(presentName);
title('ISI per region, outliers cut off');
ylabel('ISI (s)');
xlabel('Region');
ylim([-1 7]);
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(label) '.ISIboxCum.png']);

% Cell Level Burst Index
figure('Position', get(0, 'Screensize'));
boxplot(burst_boxx, burst_boxg);
xticklabels(presentName);
title('Burst Index per region');
ylabel('Burst Index');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(label) '.BurstBoxCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(perc_boxx, perc_boxg);
xticklabels(presentName);
title('Percent of units making up event per region');
ylabel('Percent of units');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(label) '.PercUntsCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(avgSpk_boxx, avgSpk_boxg);
xticklabels(presentName);
title('Average number of spikes per event per region');
ylabel('Number of spikes');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(label) '.AvgSpkCum.png']);

%
figure('Position', get(0, 'Screensize'));
boxplot(avgFR_boxx, avgFR_boxg);
xticklabels(presentName);
title('Average Firing Rate per event per region');
ylabel('Firing Rate (Hz)');
xlabel('Region');
saveas(gcf,['Z:\home\Lindsay\Barrage\cumMet\' convertStringsToChars(label) '.AvgFRCum.png']);

% cumMet.boxxISI
% cumMet.boxgISI
% 
% % cumMet.evtBin
% % cumMet.evtAvg
% 
% cumMet.boxxPerc
% cumMet.boxxgPerc
% 
% cumMet.boxxAvgSpk
% cumMet.boxgAvgSpk
% 
% cumMet.boxxAvgFR
% cumMet.boxgAvgFR
end