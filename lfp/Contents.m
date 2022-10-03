% LFP
%
% Files
%   auto_theta_cycles  - auto_theta_cycles: automatically select theta channel and calc cycles
%   bz_Filter          - Filter - Filter samples.
%   CSD                - - Calculate the 1D approximation of current source density (CSD) from LFPs
%   DetectDSpikes_v4   - DetectDSpikes_v4 - detects spikes?
%   Detrend            - Detrend - Detrend signal.
%   doICA              - [ica] = doICA(varargin)
%   eventCSD           - [ CSD ] = eventCSD (lfp, events, varargin)
%   Filter             - Filter - Filter samples.
%   FilterLFP          - FilterLFP - Filter the local field potentials, e.g. in the theta band.
%   FindDeltaWaves     - FindDeltaWaves - Find cortical delta waves (1-6Hz waves).
%   FindSpindles       - FindSpindles - Find thalamo-cortical spindles (9-17Hz oscillations).
%   FindThetaCycles    - FindThetaCycles - Find intervals that qualify as theta cycles from the lfp signal
%   gradDescCluster    - [cluass] = gradDescCluster(simMat) clusters recording sites given a pairwise 
%   Phase              - Phase - Compute instantaneous phase in signal.
%   runica             - runica() - Perform Independent Component Analysis (ICA) decomposition
%   SaveEEG            - This function will make a new lfp file (.eeg) in which the new reference
%   SineWavePeaks      - SineWavePeaks - Find peaks (or troughs) in a sine wave.
%   TransformPhaseECDF - TransformPhaseECDF - correct phases for signal nonuniformity (e.g. asymmetry)
%   WaveletSpectrogram - WaveletSpectrogram - Compute LFP wavelet spectrogram
