function behavioral_ROI_locator(varargin)
    % BEHAVIORAL_ROI_LOCATOR
    %  Locate behavioral regions of interest (ROIs) in a video for a given session.
    %
    %  BEHAVIORAL_ROI_LOCATOR 
    %  Locates behavioral ROIs in a video for a given session(s) that can integrate with the 
    % general_behavior_file. The overall goal being to label regions of interest that allow for 
    % integration within the tracking of the animal towards that ROI.
    %   Possible uses include locating objects in a maze, water spouts, or other objects that
    % the animal interacts with for a given session. The overall goal is to be able to integrate 
    % the tracking of the animal towards that ROI with physiological data.
    %
    %  NOTE: This function is dependent on the general_behavior_file and is intended to be used
    %        after the general_behavior_file has been created. Also, conversion to cm is required.
    %
    % INPUTS:
    %   [input parser]  [inputs as opiton list of values - see below]
    %   <options>       optional list of property-value pairs (see table below)
    %  =========================================================================
    %   Properties    Values
    %  -------------------------------------------------------------------------
    %  varargin:
    %  'basepath': path to the session folder. Default is pwd.
    %  'ROI_type': type of ROI to locate. Options include 'circle' or 'rectangle'. Default is circle.
    %  'ROI_number': number of ROIs to locate. Default is 1.
    %  'ROI_per_epoch': true or false. If true, will locate ROIs for each epoch. Default is false.
    %  'ROI_per_epoch_number': number of ROIs to locate per epoch. Default is 1.
    %  'ROI_save': true or false. If true, will save the ROIs. Default is false.
    %  'ROI_force': true or false. If true, will force the ROI locator to run even if the ROIs have
    %  'ROI_interaction_time': if empty, will set to false. If it exists, will set to true...
    %           and use the value as the distance to the object (in cm) to set as "interacting". Default is empty.
    %
    %
    % HLarsson 2023 - initial version
    % This program is released under the GNU General Public License (GPL) version 3.
    % It is provided free of charge and "as is" without any warranty.
    % You can redistribute it and/or modify it under the terms of the GNU General Public License.
    % 
    % TODO: 
    %
    % =========================================================================
    %
    % Parse inputs

    p = inputParser;
    addParameter(p,'ROI_type','circle',@ischar);
    addParameter(p,'ROI_number',1,@isnumeric);
    addParameter(p,'ROI_per_epoch',false,@islogical);
    addParameter(p,'ROI_per_epoch_number',1,@isnumeric);
    addParameter(p,'ROI_save',false,@islogical);
    addParameter(p,'ROI_force',false,@islogical);
    addParameter(p,'ROI_interaction_time',[],@isnumeric)

    parse(p,varargin{:});
    basepath = p.Results.basepath;
    ROI_type = p.Results.ROI_type;
    ROI_number = p.Results.ROI_number;
    ROI_per_epoch = p.Results.ROI_per_epoch;
    ROI_per_epoch_number = p.Results.ROI_per_epoch_number;
    ROI_save = p.Results.ROI_save;
    ROI_force = p.Results.ROI_force;
    ROI_interaction_time = p.Results.ROI_interaction_time;

    if ~iscell(basepaths)
        basepaths = {basepaths};
    end

    %iterate over basepaths for general behavior file, sesson file, and roi file
    for i = 1:length(basepaths)
        basepath = basepaths{i};
        basename = basenameFromBasepath(basepath);
        if exists([basepath, filesep, basename, '.animal.behavior.mat']) &&...
            ~ROI_force
            load([basepath, filesep, basename, '.animal_behavior.mat']);
        else disp(['No session file found for ', basename]);
            continue
        end
        if exists([basepath,filesep,basename,'.behavior.roi.mat']) &&...
            ~ROI_force
            load([basepath,filesep,basename,'.behavior.roi.mat']);
        else
            disp(['No roi file found for ', basename, '. Creating new roi file.']);
            continue
        end
        %load session

        