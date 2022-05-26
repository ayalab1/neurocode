% nDays = 43;
nDays = 55;

folders = cell(nDays,1);
for i=1:length(folders),
%     folders{i} = ['N:\V1test\AO52\day' num2str(i)];
    folders{i} = ['N:\V1test\V1Jean\day' num2str(i)];
end

fillMissingDatFiles = false;
fillTypes = [];
analogInputs = true;
analogChannels = 1;
digitalInputs = true;
digitalChannels = 1:4;
getAcceleration = true;
cleanArtifacts = false;
stateScore = false;
spikeSort = true;
getPos = false;
removeNoise = false;
runSummary = false;
SSD_path = 'D:\KiloSort';
done = false(nDays,1);
for i=1:nDays
    try
        basepath = folders{i}

        cd(basepath);
        f = dir('Kilosort*');
        if ~isempty(f), 
            done(i,1) = true;
            error('done');
        end

        %% Pull meta data

        % Get session names
        if strcmp(basepath(end),filesep)
            basepath = basepath(1:end-1);
        end
        [~,basename] = fileparts(basepath);

        % Get xml file from the parent folder
        [parentFolder,dayName] = fileparts(basepath);
        [~,projectName] = fileparts(parentFolder);
        parentXml = fullfile(parentFolder,[projectName '.xml']);
        xmlFile = fullfile(basepath,[dayName '.xml']);
        if ~exist(xmlFile,'file'), copyfile(parentXml,xmlFile); end
%         % Get xml file in order
%         xmlFile = checkFile('fileType','.xml','searchSubdirs',true);
%         xmlFile = xmlFile(1);
%         if ~(strcmp(xmlFile.folder,basepath)&&strcmp(xmlFile.name(1:end-4),basename))
%             copyfile([xmlFile.folder,filesep,xmlFile.name],[basepath,filesep,basename,'.xml'])
%         end

        % Check info.rhd
        % (assumes this will be the same across subsessions)
        rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
        rhdFile = rhdFile(1);
        if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
            copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
        end

        %% Make SessionInfo
        % Manually ID bad channels at this point. automating it would be good

        session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
        save([basename '.session.mat'],'session');
        %% Fill missing dat files of zeros
        if fillMissingDatFiles
            if isempty(fillTypes)
                fillTypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
            end
            for ii = 1:length(fillTypes)
                fillMissingDats('basepath',basepath,'fileType',fillTypes{ii});
            end
        end
        %% Concatenate sessions
        cd(basepath);

        disp('Concatenate session folders...');
        concatenateDats(pwd,0,0);

        %% Process additional inputs - CHECK FOR OUR LAB

        % Analog inputs
        % check the two different fucntions for delaing with analog inputs and proably rename them
        if analogInputs
            try
                if  ~isempty(analogChannels)
                    analogInp = computeAnalogInputs('analogCh',analogChannels,'saveMat',true,'fs',session.extracellular.sr);
                else
                    analogInp = computeAnalogInputs('analogCh',[],'saveMat',true,'fs',session.extracellular.sr);
                end

                % analog pulses ...
                [pulses] = getAnalogPulses('samplingRate',session.extracellular.sr);
            end
        end

        % Digital inputs
        try
            digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
        catch
            display('Error with digital inputs. This step was skipped.');
        end
        % Auxilary input
        if getAcceleration
            try
                accel = computeIntanAccel('saveMat',true);
            end
        end

        %% Make LFP

        LFPfromDat(pwd,'outFs',1250,'useGPU',false);

        % 'useGPU'=true ; gives an error if CellExplorer in the path. Need to test if
        % it is possible to remove the copy of iosr toolbox from CellExplorer

        %% Clean data  - CHECK FOR OUR LAB
        % Remove stimulation artifacts
        if cleanArtifacts && analogInputs
            [pulses] = getAnalogPulses(analogInp,'analogCh',analogChannels);
            cleanPulses(pulses.ints{1}(:));
        end

        % remove noise from data for cleaner spike sorting
        % if removeNoise
        %     NoiseRemoval(pwd); % not very well tested yet: Raly: this is bad
        % end

        %% Get brain states
        % an automatic way of flaging bad channels is needed
        if stateScore
            try
                if exist('pulses','var')
                    SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
                    thetaEpochs(pwd);
                else
                    SleepScoreMaster(pwd,'noPrompts',true); % takes lfp in base 0
                    thetaEpochs(pwd);
                end
            catch
                disp('Problem with SleepScore skyping...');
            end
        end

        %% Kilosort concatenated sessions
        if isempty(dir('KilosortGT*')) && spikeSort
            kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
            load(fullfile(kilosortFolder,'rez.mat'),'rez');
            CleanRez(rez,'savepath',kilosortFolder);
            %     PhyAutoClustering(kilosortFolder);
        end
        % move datfile to kilosort folder
        % movefile(fullfile(basepath,[basename '.dat']),fullfile(kilosortFolder,[basename '.dat']));

    catch
        warning(['Issue with session ' num2str(i)]);
    end
end
