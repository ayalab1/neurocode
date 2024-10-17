basepath = pwd;
basename = basenameFromBasepath(basepath);
load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
unit = 3;

%% ACGs
figure(1); area(cell_metrics.acg.narrow(:,unit));
title('ACG');
xlim([1 length(cell_metrics.acg.narrow(:,unit))]);
xticks([1 101 201]); xticklabels({'-50','0','50'}); xlabel('Time (ms)');
yticks([]);

%% Waveforms
wave = cell_metrics.waveforms.filt{unit};
wave_SEM = cell_metrics.waveforms.filt_std{unit};
% wave_SEM = wave_SEM/sqrt(cell_metrics.spikeCount(Pb_connection(eventIdx,1)));
figure(2);
shadedErrorBar(cell_metrics.waveforms.time{unit},wave,wave_SEM,{'LineWidth',1.5,'Color','k'},0);hold on;
title('Waveform'); xlabel('Time (ms)'); ylabel('Voltage (uV)');