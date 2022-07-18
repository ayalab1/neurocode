%% Load data
basepath=pwd;
basename = basenameFromBasepath(basepath);


load([basename,'.session.mat'])
load([basename,'.cell_metrics.cellinfo.mat'])
load([basename,'.spikes.cellinfo.mat'])
spikesCell = cell_metrics.spikes.times';
%

%%load neurons
region=["PFC"];
int_spikes= importSpikes('brainRegion', region, 'cellType', [ "Narrow Interneuron"]);
pyn_spikes= importSpikes('brainRegion', region, 'cellType', ["Pyramidal Cell"]);


%%loadpulseevents

load([basename,'.pulses.events.mat']);

hold on

PSTH_int =computePSTH(pulses,int_spikes,'duration',.3,'zscorePlot',true);
PlotHVLines(0.2,'v')

hold off
hold on
PSTH_pyn =computePSTH(pulses,pyn_spikes,'duration',.3,'zscorePlot',true);
PlotHVLines(0.2,'v')

hold off

hold on 
figure


PSTH_out_int=PSTH_int.responsecurve;
time_int=PSTH_int.time;
data_SEM_in = std(PSTH_out_int,[],2)/sqrt(size(PSTH_out_int,2)); 
    plot(time_int,mean(PSTH_out_int')','LineWidth',2); hold on; 
    plot(time_int,mean(PSTH_out_int')+data_SEM_in','--b');hold on; 
    plot(time_int,mean(PSTH_out_int')-data_SEM_in','--b');hold on;
    xline(0,'--k');hold on; ylabel('mod. index');
   
    xline(0,'--k');hold on; 

hold off
hold on
figure
PSTH_out_pyn=PSTH_pyn.responsecurve;
PSTH_out_pyn_trans=PSTH_pyn.responsecurve';

time_pyn=PSTH_pyn.time;
data_SEM_pyn = std(PSTH_out_pyn,[],2)/sqrt(size(PSTH_out_pyn,2)); 
    plot(time_pyn,mean(PSTH_out_pyn')','LineWidth',2); hold on;     
    plot(time_pyn,mean(PSTH_out_pyn')+data_SEM_pyn','--b');hold on; 
    plot(time_pyn,mean(PSTH_out_pyn')-data_SEM_pyn','--b');hold on;
    xline(0,'--k');hold on; ylabel('mod. index');
    
    xline(0,'--k');hold on; 