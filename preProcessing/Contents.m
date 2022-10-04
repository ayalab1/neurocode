% PREPROCESSING
%
% Files
%   channel_mapping            - channel_mapping: updates basename.session with brain regions per channel
%   cleanPulses                - [data] = cleanPulses(ts, varargin)
%   CleanRez                   - Remove noisy clusters from the rez file. This is an optional preprocessing
%   computeIntanAccel          - computeIntanAccel - Get accelerometer data from Intan-generated auxiliary.dat file
%   concatenateDats            - concatenateDats - Concatenate raw .dat files found in a session folder
%   concatenatedTimes          - concatenatedTimes(basepath)
%   concatenateOtherDats       - concatenateOtherDats - Function to help you concatenate auxiliary .dat files that were missed during preprocessing
%   fillMissingDats            - 
%   forceMergerOfInputFiles    - This function will produce a digitalIn.dat and analogin.dat files even if
%   getAnalogPulses            - [pul, val, dur] = getAnalogPulses(varargin)
%   hackInfo                   - patch for dealing with different session info formats - just call this to
%   helper_RemoveDatFile       - helper_RemoveDatFile- script that checks if the .dat file in the current directory matches the size of the individual subsession .dat files and removes it
%   InterpolateDat             - Replace a noisy period in a .dat file (interpolating)
%   LFPfromDat                 - perform lowpass (2 X output Fs) sinc filter on wideband data
%   multiCellMetrics           - multiCellMetrics- to visualize multiple sessions in cellExplorer GUI
%   NoiseRemoval               - Copyright (C) 2021 KM and 2022 Ralitsa Todorova
%   preprocessSession          - preprocessSession(varargin)
%   PreprocessSpikes           - PreprocessSpikes - secondary preprocessing for after manual clustering
%   ProcessBinary              - ProcessBinary - Process binary data file.
%   read_Intan_Info_Wrapper    - varargout = read_Intan_Info_Wrapper(filename, varargin)
%   read_Intan_RHD2000_file_bz - read_Intan_RHD2000_file_bz
%   ResampleBinary             - ResampleBinary - Resample binary data file.
%   script_CutDat              - move to the subfolder in question
