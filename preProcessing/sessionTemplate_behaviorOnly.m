function session = sessionTemplate_behaviorOnly(input1, varargin)
    % sessionTemplate_behaviorOnly - Create session metadata for behavior-only recordings
    %
    % This function creates a minimal session struct for recordings that only contain
    % behavioral data (digital/analog inputs, tracking) without electrophysiology data
    % (no amplifier.dat, no spike sorting needed).
    %
    % It reads metadata from:
    %   - info.rhd (Intan) or settings.xml (OpenEphys)
    %   - time.dat (to calculate duration)
    %   - digitalin.dat, analogin.dat, digitalout.dat (for channel info)
    %
    % USAGE:
    %   session = sessionTemplate_behaviorOnly(basepath)
    %   session = sessionTemplate_behaviorOnly(basepath, 'showGUI', true)
    %   session = sessionTemplate_behaviorOnly(session) % update existing session
    %
    % INPUTS:
    %   input1              - Either basepath (string) or existing session struct
    %   showGUI             - Show session GUI editor (default: false)
    %   basename            - Override basename detection
    %   addBehaviorFiles    - Auto-detect tracking/video/DLC files (default: false)
    %   nDigitalChannels    - Number of digital input channels (default: 16 for Intan)
    %   nAnalogChannels     - Number of analog input channels (default: 8 for Intan)
    %
    % OUTPUTS:
    %   session         - Session metadata structure
    %
    % EXAMPLE:
    %   session = sessionTemplate_behaviorOnly('/path/to/behavior/session');
    %   
    %   % For a system with 32 digital channels:
    %   session = sessionTemplate_behaviorOnly(basepath, 'nDigitalChannels', 32);
    %
    % See also: sessionTemplate, preprocessSession

    % Parse inputs
    p = inputParser;
    addRequired(p, 'input1', @(X) (ischar(X) && exist(X, 'dir')) || isstruct(X));
    addParameter(p, 'basename', [], @ischar);
    addParameter(p, 'showGUI', false, @islogical);
    addParameter(p, 'addBehaviorFiles', false, @islogical);
    addParameter(p, 'nDigitalChannels', 16, @isnumeric); % Intan default
    addParameter(p, 'nAnalogChannels', 8, @isnumeric);  % Intan default

    parse(p, input1, varargin{:})
    basename = p.Results.basename;
    showGUI = p.Results.showGUI;
    addBehaviorFiles = p.Results.addBehaviorFiles;
    nDigitalChannels = p.Results.nDigitalChannels;
    nAnalogChannels = p.Results.nAnalogChannels;

    % Initialize session struct and basepath
    if ischar(input1)
        basepath = input1;
        cd(basepath)
    elseif isstruct(input1)
        session = input1;

        if isfield(session.general, 'basePath') && exist(session.general.basePath, 'dir')
            basepath = session.general.basePath;
            cd(basepath)
        else
            basepath = pwd;
        end

    end

    % Get basename
    if isempty(basename)
        basename = basenameFromBasepath(basepath);
    end

    % Load existing session file if it exists
    if ~exist('session', 'var') && exist(fullfile(basepath, [basename, '.session.mat']), 'file')
        disp('Loading existing basename.session.mat file')
        session = loadSession(basepath, basename);
    elseif ~exist('session', 'var')
        session = [];
    end

    pathPieces = regexp(basepath, filesep, 'split');

    %% General metadata
    session.general.version = 5;
    session.general.name = basename;
    session.general.basePath = basepath;
    session.general.sessionType = 'Behavior';
    session.general.experimenters = ''; % Add experimenter names if known

    %% Animal metadata
    if ~isfield(session, 'animal')
        session.animal.name = pathPieces{end - 1};
        session.animal.sex = 'Unknown';
        session.animal.species = 'Unknown';
        session.animal.strain = 'Unknown';
        session.animal.geneticLine = '';
    end

    %% Detect recording system and read metadata
    rhdFile = dir(fullfile(basepath, '**', '*.rhd'));
    settingsFile = dir(fullfile(basepath, '**', 'settings.xml'));

    if ~isempty(rhdFile)
        % Intan system detected
        disp('Intan system detected (info.rhd found)')
        session = readIntanBehaviorMetadata(session, basepath, rhdFile(1));

    elseif ~isempty(settingsFile)
        % OpenEphys system detected
        disp('OpenEphys system detected (settings.xml found)')
        session = readOpenEphysBehaviorMetadata(session, basepath, settingsFile(1));

    else
        % No metadata files found - use defaults and warn user
        warning('No info.rhd or settings.xml found. Using default parameters.');
        session = setDefaultBehaviorMetadata(session, basepath);
    end

    %% Calculate duration from time.dat or data files
    session = calculateSessionDuration(session, basepath, basename);

    %% Create epochs from MergePoints or default
    session = createBehaviorEpochs(session, basepath, basename);

    %% Extract date/time from folder or files
    session = extractSessionDateTime(session, basepath, basename);

    %% Add behavioral data paths (if requested)
    if addBehaviorFiles
        session = addBehavioralPaths(session, basepath, basename);
    else
        % Remove behavioral file fields if they exist (cleanup from previous runs)
        fieldsToRemove = {'videos', 'deepLabCut', 'optitrack', 'behavioralTracking'};

        for i = 1:length(fieldsToRemove)

            if isfield(session, fieldsToRemove{i})
                session = rmfield(session, fieldsToRemove{i});
            end

        end

    end

    %% Show GUI if requested
    if showGUI
        session = gui_session(session);
    end

end

%% Helper Functions

function session = readIntanBehaviorMetadata(session, basepath, rhdFile)
    % readIntanBehaviorMetadata - Extract metadata from Intan RHD2000 files
    %
    % Reads recording parameters from an Intan info.rhd file including:
    % - Sampling rates (main and LFP)
    % - Digital input channels (board_dig_in)
    % - Analog input channels (board_adc)
    % - Auxiliary channels (accelerometer, etc.)
    % - Recording notes
    %
    % INPUTS:
    %   session     - Partially filled session struct
    %   basepath    - Path to session directory
    %   rhdFile     - Dir struct entry for the .rhd file
    %
    % OUTPUTS:
    %   session     - Session struct with Intan metadata populated
    %
    % NOTES:
    %   - Uses read_Intan_RHD2000_file_bz() to parse binary RHD format
    %   - Falls back to setDefaultBehaviorMetadata() on error
    %   - Sets extracellular.nChannels = 0 (behavior-only, no amplifier data)

    try
        % Use existing Intan reader
        rhdPath = fullfile(rhdFile.folder, rhdFile.name);

        % Read RHD file
        [amplifier_channels, notes, aux_input_channels, ~, ...
             board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels] = ...
            read_Intan_RHD2000_file_bz('basepath', rhdFile.folder);

        % Sampling rate
        session.extracellular.sr = frequency_parameters.amplifier_sample_rate;
        session.extracellular.srLfp = 1250; % Standard LFP rate (not used in behavior-only)

        % Digital channels
        if ~isempty(board_dig_in_channels)
            nDigitalChannels = length(board_dig_in_channels);
            session.inputs.nDigitalChannels = nDigitalChannels;
            session.inputs.digitalChannels = 1:nDigitalChannels;

            for i = 1:nDigitalChannels
                session.inputs.digitalChannelNames{i} = board_dig_in_channels(i).native_channel_name;
            end

        else
            session.inputs.nDigitalChannels = 0;
        end

        % Analog channels (board ADC)
        if ~isempty(board_adc_channels)
            nAnalogChannels = length(board_adc_channels);
            session.inputs.nAnalogChannels = nAnalogChannels;
            session.inputs.analogChannels = 1:nAnalogChannels;

            for i = 1:nAnalogChannels
                session.inputs.analogChannelNames{i} = board_adc_channels(i).native_channel_name;
            end

        else
            session.inputs.nAnalogChannels = 0;
        end

        % Auxiliary channels (accelerometer, etc.)
        if ~isempty(aux_input_channels)
            nAuxChannels = length(aux_input_channels);
            session.inputs.nAuxChannels = nAuxChannels;

            for i = 1:nAuxChannels
                session.inputs.auxChannelNames{i} = aux_input_channels(i).native_channel_name;
            end

        else
            session.inputs.nAuxChannels = 0;
        end

        % Notes
        if ~isempty(notes) && isstruct(notes)
            % Handle different possible field names for notes
            if isfield(notes, 'note1')
                session.general.notes = notes.note1;
            elseif isfield(notes, 'notes')
                session.general.notes = notes.notes;
            else
                % Use first available field
                noteFields = fieldnames(notes);

                if ~isempty(noteFields)
                    session.general.notes = notes.(noteFields{1});
                end

            end

        elseif ischar(notes) || isstring(notes)
            % Handle case where notes is directly a string
            session.general.notes = char(notes);
        end

        % Set minimal extracellular fields (even though no spikes)
        session.extracellular.nChannels = 0; % No amplifier channels in behavior-only
        session.extracellular.nElectrodeGroups = 0;
        session.extracellular.nSpikeGroups = 0;
        session.extracellular.precision = 'int16'; % Standard precision for raw data
        session.extracellular.leastSignificantBit = 0.195; % Intan default (in µV)

        disp(['  Sampling rate: ', num2str(session.extracellular.sr), ' Hz'])
        disp(['  Digital channels: ', num2str(session.inputs.nDigitalChannels)])
        disp(['  Analog channels: ', num2str(session.inputs.nAnalogChannels)])

    catch ME
        warning(['Failed to read Intan metadata: ', ME.message]);
        session = setDefaultBehaviorMetadata(session, basepath);
    end

end

function session = readOpenEphysBehaviorMetadata(session, basepath, settingsFile)
    % readOpenEphysBehaviorMetadata - Extract metadata from OpenEphys settings.xml
    %
    % Parses OpenEphys XML configuration to extract recording parameters.
    % Attempts to find sampling rate from SIGNALCHAIN->PROCESSOR->EDITOR hierarchy.
    % Falls back to detectChannelsFromFiles() for channel detection since XML
    % may not contain complete channel information.
    %
    % INPUTS:
    %   session         - Partially filled session struct
    %   basepath        - Path to session directory
    %   settingsFile    - Dir struct entry for settings.xml
    %
    % OUTPUTS:
    %   session         - Session struct with OpenEphys metadata populated
    %
    % NOTES:
    %   - Default sampling rate is 30000 Hz if not found in XML
    %   - Uses xml2struct() for XML parsing
    %   - Calls detectChannelsFromFiles() to count actual data channels
    %   - Falls back to setDefaultBehaviorMetadata() on error

    try
        % Parse XML file
        settingsPath = fullfile(settingsFile.folder, settingsFile.name);
        xmlStruct = xml2struct(settingsPath);

        % Extract sampling rate (typically in SIGNALCHAIN > PROCESSOR > EDITOR)
        % OpenEphys default is 30000 Hz
        session.extracellular.sr = 30000; % Default, will try to parse from XML

        % Try to extract actual sampling rate from XML
        try

            if isfield(xmlStruct, 'SETTINGS') && isfield(xmlStruct.SETTINGS, 'SIGNALCHAIN')
                signalChain = xmlStruct.SETTINGS.SIGNALCHAIN;

                if isfield(signalChain, 'PROCESSOR')
                    processors = signalChain.PROCESSOR;

                    if ~iscell(processors)
                        processors = {processors};
                    end

                    for i = 1:length(processors)

                        if isfield(processors{i}, 'Attributes') && isfield(processors{i}.Attributes, 'name')

                            if contains(processors{i}.Attributes.name, 'Sources', 'IgnoreCase', true)

                                if isfield(processors{i}, 'EDITOR') && isfield(processors{i}.EDITOR, 'Attributes')

                                    if isfield(processors{i}.EDITOR.Attributes, 'SampleRate')
                                        session.extracellular.sr = str2double(processors{i}.EDITOR.Attributes.SampleRate);
                                    end

                                end

                            end

                        end

                    end

                end

            end

        catch
            warning('Could not parse sampling rate from settings.xml, using default 30000 Hz');
        end

        session.extracellular.srLfp = 1250;

        % Set minimal extracellular fields
        session.extracellular.nChannels = 0;
        session.extracellular.nElectrodeGroups = 0;
        session.extracellular.nSpikeGroups = 0;
        session.extracellular.precision = 'int16'; % Standard precision for raw data
        session.extracellular.leastSignificantBit = 0.195; % Default (Intan = 0.195, Amplipex = 0.3815)

        % Detect digital/analog channels from actual files
        session = detectChannelsFromFiles(session, basepath, nDigitalChannels, nAnalogChannels);

        disp(['  Sampling rate: ', num2str(session.extracellular.sr), ' Hz'])
        disp(['  Digital channels detected: ', num2str(session.inputs.nDigitalChannels)])
        disp(['  Analog channels detected: ', num2str(session.inputs.nAnalogChannels)])

    catch ME
        warning(['Failed to read OpenEphys metadata: ', ME.message]);
        session = setDefaultBehaviorMetadata(session, basepath);
    end

end

function session = setDefaultBehaviorMetadata(session, basepath)
    % setDefaultBehaviorMetadata - Apply fallback metadata when no info files exist
    %
    % Used when neither info.rhd (Intan) nor settings.xml (OpenEphys) are found.
    % Sets conservative default values and relies on detectChannelsFromFiles()
    % to determine actual channel counts from data files.
    %
    % INPUTS:
    %   session     - Partially filled session struct
    %   basepath    - Path to session directory
    %
    % OUTPUTS:
    %   session     - Session struct with default metadata
    %
    % DEFAULTS:
    %   - Sampling rate: 20000 Hz (common Intan default)
    %   - LFP rate: 1250 Hz
    %   - No electrode groups or spike groups (behavior-only)
    %   - Channel counts determined from actual files
    %
    % NOTES:
    %   - Displays warning that user should verify parameters
    %   - Should only be used as last resort

    warning('Using default behavior metadata. Please verify sampling rate and channel counts!');

    session.extracellular.sr = 20000; % Intan default
    session.extracellular.srLfp = 1250;
    session.extracellular.nChannels = 0;
    session.extracellular.nElectrodeGroups = 0;
    session.extracellular.nSpikeGroups = 0;

    % Detect channels from files
    session = detectChannelsFromFiles(session, basepath, nDigitalChannels, nAnalogChannels);

end

function session = detectChannelsFromFiles(session, basepath, defaultDigitalChannels, defaultAnalogChannels)
    % detectChannelsFromFiles - Determine channel counts from actual data files
    %
    % Searches for digitalin.dat, analogin.dat, and digitalout.dat files to
    % determine what input channels are present. Uses configurable defaults
    % for channel counts since files don't self-describe their channel count.
    %
    % INPUTS:
    %   session                  - Session struct being built
    %   basepath                 - Path to session directory
    %   defaultDigitalChannels   - Number of digital input channels (default: 16 for Intan)
    %   defaultAnalogChannels    - Number of analog input channels (default: 8 for Intan)
    %
    % OUTPUTS:
    %   session                  - Session struct with populated channel information:
    %                              - inputs.nDigitalChannels, digitalChannels, digitalChannelNames
    %                              - inputs.nAnalogChannels, analogChannels, analogChannelNames
    %                              - inputs.hasDigitalOut
    %
    % FILE FORMATS:
    %   digitalin.dat  - uint16, multiplexed digital channels
    %   analogin.dat   - uint16, multiplexed analog channels
    %   digitalout.dat - uint16, digital output channels
    %
    % NOTES:
    %   - Searches recursively with dir(fullfile(basepath, '**', filename))
    %   - Channel names are auto-generated (e.g., 'Digital-0', 'Analog-0')
    %   - Defaults can be overridden via function parameters
    
    % Check for digitalin.dat
    digitalFiles = dir(fullfile(basepath, '**', 'digitalin.dat'));

    if ~isempty(digitalFiles)
        digitalFile = digitalFiles(1);
        fileInfo = dir(fullfile(digitalFile.folder, digitalFile.name));
        % Use the configurable default
        session.inputs.nDigitalChannels = defaultDigitalChannels;
        session.inputs.digitalChannels = 1:defaultDigitalChannels;
        
        % Generate default channel names
        for i = 1:defaultDigitalChannels
            session.inputs.digitalChannelNames{i} = sprintf('Digital-%d', i-1);
        end
    else
        session.inputs.nDigitalChannels = 0;
    end

    % Check for analogin.dat
    analogFiles = dir(fullfile(basepath, '**', 'analogin.dat'));

    if ~isempty(analogFiles)
        analogFile = analogFiles(1);
        fileInfo = dir(fullfile(analogFile.folder, analogFile.name));
        % Use the configurable default
        session.inputs.nAnalogChannels = defaultAnalogChannels;
        session.inputs.analogChannels = 1:defaultAnalogChannels;
        
        % Generate default channel names
        for i = 1:defaultAnalogChannels
            session.inputs.analogChannelNames{i} = sprintf('Analog-%d', i-1);
        end
    else
        session.inputs.nAnalogChannels = 0;
    end

    % Check for digitalout.dat
    digitalOutFiles = dir(fullfile(basepath, '**', 'digitalout.dat'));

    if ~isempty(digitalOutFiles)
        session.inputs.hasDigitalOut = true;
    else
        session.inputs.hasDigitalOut = false;
    end

end

function session = calculateSessionDuration(session, basepath, basename)
    % calculateSessionDuration - Determine recording duration from data files
    %
    % Attempts to calculate total session duration using available files,
    % prioritizing time.dat for accuracy, then falling back to data files.
    %
    % INPUTS:
    %   session     - Session struct being built
    %   basepath    - Path to session directory
    %   basename    - Base name for session files (unused but kept for compatibility)
    %
    % OUTPUTS:
    %   session     - Session struct with populated:
    %                 - extracellular.nSamples (total samples)
    %                 - general.duration (seconds)
    %
    % FILE PRIORITY:
    %   1. time.dat       - Most accurate (contains actual timestamps)
    %   2. digitalin.dat  - Good fallback (file size / 2 bytes = samples)
    %   3. analogin.dat   - Last resort (file size / (nChannels * 2 bytes) = samples)
    %
    % NOTES:
    %   - time.dat contains int32 timestamps in sample units
    %   - Duration = nSamples / samplingRate
    %   - Searches recursively for files in subdirectories
    %   - Displays warning if no usable files found

    % Try to use time.dat first
    timeFiles = dir(fullfile(basepath, '**', 'time.dat'));

    if ~isempty(timeFiles)
        timeFile = fullfile(timeFiles(1).folder, timeFiles(1).name);

        try
            fid = fopen(timeFile, 'r');
            timeStamps = fread(fid, inf, 'int32');
            fclose(fid);

            if ~isempty(timeStamps)
                % time.dat contains timestamps in samples
                session.extracellular.nSamples = double(max(timeStamps));
                session.general.duration = session.extracellular.nSamples / session.extracellular.sr;
                disp(['  Duration (from time.dat): ', num2str(session.general.duration), ' seconds'])
                return
            end

        catch
            warning('Could not read time.dat');
        end

    end

    % Fall back to digitalin.dat or analogin.dat
    digitalFiles = dir(fullfile(basepath, '**', 'digitalin.dat'));

    if ~isempty(digitalFiles)
        digitalFile = fullfile(digitalFiles(1).folder, digitalFiles(1).name);
        fileInfo = dir(digitalFile);
        % digitalin.dat: uint16, 2 bytes per sample
        nSamples = fileInfo.bytes / 2;
        session.extracellular.nSamples = nSamples;
        session.general.duration = nSamples / session.extracellular.sr;
        disp(['  Duration (from digitalin.dat): ', num2str(session.general.duration), ' seconds'])
        return
    end

    % Try analogin.dat
    analogFiles = dir(fullfile(basepath, '**', 'analogin.dat'));

    if ~isempty(analogFiles)
        analogFile = fullfile(analogFiles(1).folder, analogFiles(1).name);
        fileInfo = dir(analogFile);
        % analogin.dat: uint16, 8 channels, 2 bytes per sample per channel
        nSamples = fileInfo.bytes / (session.inputs.nAnalogChannels * 2);
        session.extracellular.nSamples = nSamples;
        session.general.duration = nSamples / session.extracellular.sr;
        disp(['  Duration (from analogin.dat): ', num2str(session.general.duration), ' seconds'])
        return
    end

    % If nothing found
    warning('Could not calculate session duration. No time.dat, digitalin.dat, or analogin.dat found.');
    session.extracellular.nSamples = 0;
    session.general.duration = 0;

end

function s = xml2struct(xmlFile)
    % xml2struct - Convert XML file to MATLAB struct
    %
    % Recursively parses XML document into nested struct with element names
    % as field names and attributes stored in .Attributes subfields.
    %
    % INPUTS:
    %   xmlFile     - Path to XML file
    %
    % OUTPUTS:
    %   s           - Struct representation of XML hierarchy
    %
    % EXAMPLE:
    %   For XML: <SETTINGS><PROCESSOR name="test" rate="30000">...</PROCESSOR></SETTINGS>
    %   Output:  s.SETTINGS.PROCESSOR.Attributes.name = 'test'
    %            s.SETTINGS.PROCESSOR.Attributes.rate = '30000'
    %
    % NOTES:
    %   - Uses MATLAB's built-in xmlread() for parsing
    %   - Calls parseChildNodes() recursively for nested elements
    %   - Replaces '-' and ':' in element names with '_'
    %   - Element attributes stored in .Attributes substruct
    %   - Basic implementation; consider external tools for complex XML

    try
        tree = xmlread(xmlFile);
        s = parseChildNodes(tree);
    catch ME
        error(['Failed to parse XML file: ', ME.message]);
    end

end

function children = parseChildNodes(theNode)
    % parseChildNodes - Recursively parse XML DOM nodes into struct
    %
    % Helper function for xml2struct(). Traverses XML document tree and
    % converts nodes and their attributes into nested MATLAB struct.
    %
    % INPUTS:
    %   theNode     - XML DOM node object (from xmlread)
    %
    % OUTPUTS:
    %   children    - Struct containing parsed child elements and attributes
    %
    % PARSING RULES:
    %   - Only ELEMENT_NODE types are processed (ignores text, comments, etc.)
    %   - Node names with '-' or ':' are converted to '_' for valid fieldnames
    %   - Node attributes stored in .Attributes substruct
    %   - Recursively processes nested child nodes
    %
    % NOTES:
    %   - Called internally by xml2struct()
    %   - Uses Java DOM API (getChildNodes, getAttributes, etc.)

    children = struct;

    if theNode.hasChildNodes
        childNodes = theNode.getChildNodes;
        numChildNodes = childNodes.getLength;

        for count = 1:numChildNodes
            theChild = childNodes.item(count - 1);

            if theChild.getNodeType == theChild.ELEMENT_NODE
                childName = char(theChild.getNodeName);
                childName = strrep(childName, '-', '_');
                childName = strrep(childName, ':', '_');

                % Get attributes
                if theChild.hasAttributes
                    attributes = theChild.getAttributes;
                    numAttributes = attributes.getLength;
                    attributeStruct = struct;

                    for i = 1:numAttributes
                        attr = attributes.item(i - 1);
                        attrName = char(attr.getName);
                        attrValue = char(attr.getValue);
                        attributeStruct.(attrName) = attrValue;
                    end

                    children.(childName).Attributes = attributeStruct;
                end

                % Recurse on child nodes
                if theChild.hasChildNodes
                    childData = parseChildNodes(theChild);

                    if ~isempty(fieldnames(childData))
                        children.(childName) = childData;
                    end

                end

            end

        end

    end

end

function session = createBehaviorEpochs(session, basepath, basename)
    % createBehaviorEpochs - Define temporal epochs for the session
    %
    % Creates epoch definitions from available information, prioritizing
    % MergePoints.events.mat (from concatenateDats), then subsession folders,
    % finally falling back to single default epoch.
    %
    % INPUTS:
    %   session     - Partially filled session struct
    %   basepath    - Path to session directory
    %   basename    - Base name for session files
    %
    % OUTPUTS:
    %   session     - Session struct with populated epochs cell array
    %                 Each epoch contains: name, startTime, stopTime,
    %                 behavioralParadigm, environment, manipulation
    %
    % EPOCH SOURCES (in priority order):
    %   1. MergePoints.events.mat - Most accurate (from concatenateDats)
    %   2. Subsession folders (session_01, session_02, etc.) - Approximate
    %   3. Single default epoch - Last resort
    %
    % BEHAVIORAL INFERENCE:
    %   - Attempts to infer behavioralParadigm from folder names
    %   - Keywords: 'sleep'/'rest' → Sleep, 'maze'/'task' → Maze, 'open'/'field' → Open field
    %
    % NOTES:
    %   - Skips update if existing epochs already match MergePoints
    %   - Warns user if MergePoints missing with multiple subsessions
    %   - Subsession epoch times are approximate until concatenateDats runs

    % Check for MergePoints file (created by concatenateDats)
    mergePointsFile = fullfile(basepath, [basename, '.MergePoints.events.mat']);

    if exist(mergePointsFile, 'file')

        try
            load(mergePointsFile, 'MergePoints');

            % Clear existing epochs if they don't match MergePoints
            shouldUpdate = true;

            if isfield(session, 'epochs') && ~isempty(session.epochs)
                % Check if epochs already match MergePoints
                if length(session.epochs) == size(MergePoints.foldernames, 2)
                    % Check if timestamps match
                    allMatch = true;

                    for i = 1:length(session.epochs)

                        if session.epochs{i}.startTime ~= MergePoints.timestamps(i, 1) || ...
                                session.epochs{i}.stopTime ~= MergePoints.timestamps(i, 2)
                            allMatch = false;
                            break;
                        end

                    end

                    if allMatch
                        disp('  Epochs already match MergePoints, skipping update');
                        shouldUpdate = false;
                    end

                end

            end

            if shouldUpdate
                disp('Creating epochs from MergePoints...');
                % Clear old epochs
                session.epochs = {};

                for i = 1:size(MergePoints.foldernames, 2)
                    session.epochs{i}.name = MergePoints.foldernames{i};
                    session.epochs{i}.startTime = MergePoints.timestamps(i, 1);
                    session.epochs{i}.stopTime = MergePoints.timestamps(i, 2);

                    % Add behavioral context if available
                    session.epochs{i}.behavioralParadigm = 'Unknown';
                    session.epochs{i}.environment = 'Unknown';
                    session.epochs{i}.manipulation = 'None';

                    % Try to infer from folder name
                    folderName = lower(MergePoints.foldernames{i});

                    if contains(folderName, 'sleep') || contains(folderName, 'rest')
                        session.epochs{i}.behavioralParadigm = 'Sleep';
                        session.epochs{i}.environment = 'Sleep box';
                    elseif contains(folderName, 'maze') || contains(folderName, 'task') || contains(folderName, 'session')
                        session.epochs{i}.behavioralParadigm = 'Maze';
                        session.epochs{i}.environment = 'Maze';
                    elseif contains(folderName, 'open') || contains(folderName, 'field')
                        session.epochs{i}.behavioralParadigm = 'Open field';
                        session.epochs{i}.environment = 'Open field';
                    end
                end

                disp(['  Created ', num2str(length(session.epochs)), ' epochs from MergePoints']);
            end

            return

        catch ME
            warning(['Could not load MergePoints: ', ME.message]);
        end

    end

    % If no MergePoints, try to detect subsession folders
    subsessionFolders = dir(fullfile(basepath, 'session*'));
    subsessionFolders = subsessionFolders([subsessionFolders.isdir]);

    if length(subsessionFolders) >= 2
        % Multiple subsessions found - create epochs from them
        disp(['Found ', num2str(length(subsessionFolders)), ' subsession folders']);
        disp('WARNING: MergePoints.events.mat not found. Please run concatenateDats first!');
        disp('Creating placeholder epochs from subsession folders...');

        session.epochs = {};
        cumulativeTime = 0;

        for i = 1:length(subsessionFolders)
            session.epochs{i}.name = subsessionFolders(i).name;
            session.epochs{i}.startTime = cumulativeTime;

            % Try to get duration from digitalin.dat in subsession folder
            subDatFile = fullfile(basepath, subsessionFolders(i).name, 'digitalin.dat');

            if exist(subDatFile, 'file')
                fileInfo = dir(subDatFile);
                % digitalin.dat: uint16, 2 bytes per sample
                nSamples = fileInfo.bytes / 2;
                duration = nSamples / session.extracellular.sr;
                session.epochs{i}.stopTime = cumulativeTime + duration;
                cumulativeTime = cumulativeTime + duration;
            else
                session.epochs{i}.stopTime = inf;
            end

            % Add behavioral context
            session.epochs{i}.behavioralParadigm = 'Unknown';
            session.epochs{i}.environment = 'Unknown';
            session.epochs{i}.manipulation = 'None';

            % Infer from folder name
            folderName = lower(subsessionFolders(i).name);

            if contains(folderName, 'sleep') || contains(folderName, 'rest')
                session.epochs{i}.behavioralParadigm = 'Sleep';
                session.epochs{i}.environment = 'Sleep box';
            elseif contains(folderName, 'maze') || contains(folderName, 'task') || contains(folderName, 'session')
                session.epochs{i}.behavioralParadigm = 'Maze';
                session.epochs{i}.environment = 'Maze';
            elseif contains(folderName, 'open') || contains(folderName, 'field')
                session.epochs{i}.behavioralParadigm = 'Open field';
                session.epochs{i}.environment = 'Open field';
            end

        end

        disp('  NOTE: Epoch times are approximate. Run concatenateDats for accurate timestamps.');
        return
    end

    % No MergePoints and no subsessions - check if we need a default epoch
    if ~isfield(session, 'epochs') || isempty(session.epochs)
        session = createDefaultEpoch(session, basename);
    else
        disp('  Using existing epochs (no MergePoints or subsessions found)');
    end

end

function session = createDefaultEpoch(session, basename)
    % createDefaultEpoch - Create single epoch spanning entire session
    %
    % Used when no MergePoints or subsession folders exist. Creates one
    % epoch covering the full recording duration with placeholder metadata.
    %
    % INPUTS:
    %   session     - Partially filled session struct
    %   basename    - Base name for session (used as epoch name)
    %
    % OUTPUTS:
    %   session     - Session struct with single epoch:
    %                 - name: basename
    %                 - startTime: 0
    %                 - stopTime: session.general.duration (or inf if unknown)
    %                 - behavioralParadigm: 'Unknown'
    %                 - environment: 'Unknown'
    %                 - manipulation: 'None'
    %
    % NOTES:
    %   - Used as fallback when no epoch information available
    %   - stopTime set to inf if duration not yet calculated

    session.epochs{1}.name = basename;
    session.epochs{1}.startTime = 0;

    if isfield(session.general, 'duration') && session.general.duration > 0
        session.epochs{1}.stopTime = session.general.duration;
    else
        session.epochs{1}.stopTime = inf; % Will be updated when duration is known
    end

    session.epochs{1}.behavioralParadigm = 'Unknown';
    session.epochs{1}.environment = 'Unknown';
    session.epochs{1}.manipulation = 'None';

    disp('  Created default epoch for entire session');

end

function session = extractSessionDateTime(session, basepath, ~)
    % extractSessionDateTime - Determine recording date and time
    %
    % Attempts to extract recording timestamp from folder naming conventions
    % or file metadata. Supports common Intan and OpenEphys naming patterns.
    %
    % INPUTS:
    %   session     - Partially filled session struct
    %   basepath    - Path to session directory
    %   (unused)    - Placeholder for basename (not needed)
    %
    % OUTPUTS:
    %   session     - Session struct with populated:
    %                 - general.date (format: 'YYYY-MM-DD')
    %                 - general.time (format: 'HH:MM:SS')
    %
    % SUPPORTED FORMATS:
    %   Intan:      basename_YYMMDD_HHMMSS (e.g., session_210315_143022)
    %   OpenEphys:  basename_YYYY-MM-DD_HH-MM-SS (e.g., session_2021-03-15_14-30-22)
    %
    % FALLBACK ORDER:
    %   1. Folder name pattern matching
    %   2. .rhd file modification timestamp
    %   3. Current date (last resort)
    %
    % NOTES:
    %   - Uses regex to parse common date/time patterns
    %   - Assumes YY dates are 2000+
    %   - Warns when falling back to current date

    % Try to extract from folder name first (common formats)
    pathParts = strsplit(basepath, filesep);
    sessionFolder = pathParts{end};

    % Common date/time patterns
    % Intan: basename_YYMMDD_HHMMSS
    % OpenEphys: basename_YYYY-MM-DD_HH-MM-SS
    % General: basename_date

    % Try Intan format: _YYMMDD_HHMMSS
    intanPattern = '_(\d{6})_(\d{6})';
    tokens = regexp(sessionFolder, intanPattern, 'tokens');

    if ~isempty(tokens)
        dateStr = tokens{1}{1}; % YYMMDD
        timeStr = tokens{1}{2}; % HHMMSS

        % Parse date
        year = 2000 + str2double(dateStr(1:2));
        month = str2double(dateStr(3:4));
        day = str2double(dateStr(5:6));

        % Parse time
        hour = str2double(timeStr(1:2));
        minute = str2double(timeStr(3:4));
        second = str2double(timeStr(5:6));

        session.general.date = sprintf('%04d-%02d-%02d', year, month, day);
        session.general.time = sprintf('%02d:%02d:%02d', hour, minute, second);

        disp(['  Extracted date: ', session.general.date, ' ', session.general.time]);
        return
    end

    % Try OpenEphys format: _YYYY-MM-DD_HH-MM-SS
    oePattern = '_(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})';
    tokens = regexp(sessionFolder, oePattern, 'tokens');

    if ~isempty(tokens)
        year = str2double(tokens{1}{1});
        month = str2double(tokens{1}{2});
        day = str2double(tokens{1}{3});
        hour = str2double(tokens{1}{4});
        minute = str2double(tokens{1}{5});
        second = str2double(tokens{1}{6});

        session.general.date = sprintf('%04d-%02d-%02d', year, month, day);
        session.general.time = sprintf('%02d:%02d:%02d', hour, minute, second);

        disp(['  Extracted date: ', session.general.date, ' ', session.general.time]);
        return
    end

    % Fall back to file modification time
    rhdFiles = dir(fullfile(basepath, '**', '*.rhd'));

    if ~isempty(rhdFiles)
        fileDate = datetime(rhdFiles(1).date, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        session.general.date = datestr(fileDate, 'yyyy-mm-dd');
        session.general.time = datestr(fileDate, 'HH:MM:SS');
        disp(['  Using file timestamp: ', session.general.date, ' ', session.general.time]);
        return
    end

    % Last resort: use current date
    session.general.date = datestr(now, 'yyyy-mm-dd');
    disp(['  Could not extract date from folder or files. Using current date: ', session.general.date]);

end

function session = addBehavioralPaths(session, basepath, basename)
    % addBehavioralPaths - Detect and catalog behavioral data files
    %
    % Searches for tracking, video, DeepLabCut, and Optitrack files and
    % adds their paths to the session struct. Useful for automatic discovery
    % of processed behavioral data.
    %
    % INPUTS:
    %   session     - Partially filled session struct
    %   basepath    - Path to session directory
    %   basename    - Base name for session files
    %
    % OUTPUTS:
    %   session     - Session struct with populated fields:
    %                 - behavioralTracking: Cell array of tracking file names
    %                 - videos: Cell array of video file paths (relative)
    %                 - deepLabCut: Cell array of DLC CSV file paths
    %                 - optitrack: Cell array of Optitrack CSV file paths
    %
    % DETECTED FILES:
    %   Tracking:   *.position.behavior.mat, *.tracking.behavior.mat, *.animal.behavior.mat
    %   Videos:     *.mp4, *.avi, *.mov, *.mkv (recursive search)
    %   DeepLabCut: *DLC*.csv (recursive search)
    %   Optitrack:  *optitrack*.csv (recursive search)
    %
    % NOTES:
    %   - Stores relative paths (from basepath) for portability
    %   - Searches recursively in subdirectories
    %   - Only called if 'addBehaviorFiles' parameter is true
    %   - Displays summary of found files

    % Check for tracking files
    trackingFiles = {
                     [basename, '.position.behavior.mat']
                     [basename, '.tracking.behavior.mat']
                     [basename, '.animal.behavior.mat']
                     };

    for i = 1:length(trackingFiles)

        if exist(fullfile(basepath, trackingFiles{i}), 'file')

            if ~isfield(session, 'behavioralTracking')
                session.behavioralTracking = {};
            end

            [~, name, ~] = fileparts(trackingFiles{i});
            session.behavioralTracking{end + 1} = name;
            disp(['  Found tracking file: ', trackingFiles{i}]);
        end

    end

    % Check for video files
    videoExtensions = {'*.mp4', '*.avi', '*.mov', '*.mkv'};

    for i = 1:length(videoExtensions)
        videoFiles = dir(fullfile(basepath, '**', videoExtensions{i}));

        if ~isempty(videoFiles)

            if ~isfield(session, 'videos')
                session.videos = {};
            end

            for j = 1:length(videoFiles)
                videoPath = fullfile(videoFiles(j).folder, videoFiles(j).name);
                % Store relative path
                session.videos{end + 1} = strrep(videoPath, [basepath, filesep], '');
            end

            disp(['  Found ', num2str(length(videoFiles)), ' video files (', videoExtensions{i}, ')']);
        end

    end

    % Check for DeepLabCut output
    dlcFiles = dir(fullfile(basepath, '**', '*DLC*.csv'));

    if ~isempty(dlcFiles)

        if ~isfield(session, 'deepLabCut')
            session.deepLabCut = {};
        end

        for i = 1:length(dlcFiles)
            dlcPath = fullfile(dlcFiles(i).folder, dlcFiles(i).name);
            session.deepLabCut{end + 1} = strrep(dlcPath, [basepath, filesep], '');
        end

        disp(['  Found ', num2str(length(dlcFiles)), ' DeepLabCut files']);
    end

    % Check for Optitrack CSV
    optitrackFiles = dir(fullfile(basepath, '**', '*optitrack*.csv'));

    if ~isempty(optitrackFiles)

        if ~isfield(session, 'optitrack')
            session.optitrack = {};
        end

        for i = 1:length(optitrackFiles)
            optiPath = fullfile(optitrackFiles(i).folder, optitrackFiles(i).name);
            session.optitrack{end + 1} = strrep(optiPath, [basepath, filesep], '');
        end

        disp(['  Found ', num2str(length(optitrackFiles)), ' Optitrack files']);
    end

end
