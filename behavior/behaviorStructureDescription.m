%% AYA lab standard behavior structure 

% Here is a description the of the basic behavior structure(animal.behavior.mat)
% This should be flexible enought to accomodate all tasks that we do in the
% lab. In addition to the standard fields you can also add custom ones. 
% It is based on the CellExplorer data format 
% (https://cellexplorer.org/datastructure/data-structure-and-format/#behavior)

%% The pipeline is the following:
% 1- Generate postion tracking using you favourite method (e.g. LED tracking, DeepLabCut, etc.)
% 2- Run general_behavior_file.m to generate the basic animal.beahvior.mat structure
% 3- Run the specific function for your particular task (e.g.
%   linearTrackBehavior.mat) to further populate the behavior.mat structure
%   with your task specific information 

%% Behavior structure 
% A MATLAB struct behaviorName stored in a .mat file: basename.animal.behavior.mat
% with the following fields:

%   timestamps:     array of timestamps that match the data subfields (in seconds).
%   sr:             sampling rate (Hz).
%   processinginfo: a struct with information about how the .mat file was generated 
%                   including: function name of the function, version, date, parameters.

%   position:       .x, .y, .z spatial position and .linearized (a projection of 2D postions
%                   into a 1 dimensional representation). Default units: cm.
%   speed:          a 1D representation of the running speed (cm/s).
%   acceleration:   a 1D representation of the acceleration (cm^2/s).

%   epochs:         cell array for each sub-session (e.g. sleep, task1, task2, etc).
%                   This will be typically generated from .MergePoints.mat.
%                   Contains:.startTime,.stopTime, and strings .name,.enviroment,.behaviroalParadigm
%   trials:         [start stop] for all behaviroal trials (e.g. laps in a linear maze)
%   trialID:        nTrials x 1 (or more). Vector (or matrix) with a number
%                   for each trial denoting the type of trial (e.g. inbound, outboud, etc.)
%   trialIDnames:   string array with the names that correspond to the different 
%                   numbers in .trialID
%   states:         [start stop]. This is a way to further subdivide trials
%                   or behavioral epochs into things like interaction times, regions of the
%                   maze (e.g. central/side arms), etc. 
%   stateID:        nstates x 1 (or more). Vector (or matrix) with a number
%                   for each state denoting the type of state 
%   stateIDnames:   string array with the names that correspond to the different 
%                   numbers in .stateID
