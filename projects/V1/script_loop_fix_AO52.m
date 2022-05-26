
toPreprocess = {%'N:\V1test\AO52\day14' % no data
%     'N:\V1test\AO52\day44'
%     'N:\V1test\AO52\day45'
%     'N:\V1test\AO52\day46'
%     'N:\V1test\AO52\day47'
%     'N:\V1test\AO52\day48'
%     'N:\V1test\AO52\day49'
%     'N:\V1test\AO52\day50'
%     'N:\V1test\AO52\day51'
%     'N:\V1test\AO52\day52'
%     'N:\V1test\AO52\day53'
%     'N:\V1test\AO52\day54'
    };

fillMissingDatFiles = true;
fillTypes = [];
analogInputs = true;
analogChannels = 1;
digitalInputs = true;
digitalChannels = 1:4;
getAcceleration = true;
cleanArtifacts = false;
stateScore = true;
spikeSort = true;
getPos = false;
removeNoise = false;
runSummary = false;
SSD_path = 'D:\KiloSort';
done = false(nDays,1);

%% First step: preprocess (merging dats) sessions that weren't done before:
for i=1:length(toPreprocess)
    basepath = toPreprocess{i}
    [~,basename] = fileparts(basepath);

    cd(basepath);

    disp([datestr(clock) ': Launching preprossing for ' basename]);
    % Pull meta data
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
    rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
    rhdFile = rhdFile(1);
    if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
        copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
    end

    session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
    save([basename '.session.mat'],'session');
    % Fill missing dat files of zeros
    if fillMissingDatFiles
        if isempty(fillTypes)
            fillTypes = {'analogin';'digitalin';'auxiliary';'time'};
        end
        for ii = 1:length(fillTypes)
            fillMissingDats('basepath',basepath,'fileType',fillTypes{ii});
        end
    end
    % Concatenate sessions
    cd(basepath);

    disp('Concatenate session folders...');
    concatenateDats(pwd,0,0);

    % Process additional inputs - CHECK FOR OUR LAB
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

    % Make LFP
    LFPfromDat(pwd,'outFs',1250,'useGPU',false);

    % Remove stimulation artifacts
    if cleanArtifacts && analogInputs
        [pulses] = getAnalogPulses(analogInp,'analogCh',analogChannels);
        cleanPulses(pulses.ints{1}(:));
    end

    % later steps removed
end

%% Second step, clean dat files of sessions where this was not done

toClean = {%'N:\V1test\AO52\day1'
%     'N:\V1test\AO52\day10'
%     'N:\V1test\AO52\day12'
%     'N:\V1test\AO52\day13'
%     'N:\V1test\AO52\day15'
%     'N:\V1test\AO52\day16'
%     'N:\V1test\AO52\day17'
%     'N:\V1test\AO52\day18'
%     'N:\V1test\AO52\day19'
%     'N:\V1test\AO52\day20'
%     'N:\V1test\AO52\day21'
%     'N:\V1test\AO52\day22'
%     'N:\V1test\AO52\day24'
%     'N:\V1test\AO52\day25'
%     'N:\V1test\AO52\day26'
%     'N:\V1test\AO52\day38'
%     'N:\V1test\AO52\day39'
%     'N:\V1test\AO52\day40'
%     'N:\V1test\AO52\day41'
%     'N:\V1test\AO52\day42'
%     'N:\V1test\AO52\day43'
%     'N:\V1test\AO52\day44'
%     'N:\V1test\AO52\day45'
%     'N:\V1test\AO52\day46'
%     'N:\V1test\AO52\day47'
%     'N:\V1test\AO52\day48'
%     'N:\V1test\AO52\day49'
%     'N:\V1test\AO52\day50'
%     'N:\V1test\AO52\day51'
%     'N:\V1test\AO52\day52'
%     'N:\V1test\AO52\day53'
%     'N:\V1test\AO52\day54'
    };

for i=1:length(toClean)
    basepath = toClean{i};
    [~,basename] = fileparts(basepath);
    disp([datestr(clock) ': Cleaning dat-file for ' basename]);
    try
    % Apply AO52-specific clean:
    rejectChannels = 1+[8 12 14 15]; % for AO52
    nChannels = 64;

    basename = basenameFromBasepath(basepath);
    datFile = [basepath,filesep, basename, '.dat'];
    m = memmapfile(datFile, 'Format','int16','Writable',true);
    data = reshape(m.data,64,[]);
    nSamples = size(data,2);
    okChannels = find(~ismember((1:size(data,1))',rejectChannels));
    okChannels = okChannels(1:2:end); % take every other channel, no need to saturate memory
    signal = mean(data(okChannels,:))';
    bad = [false; abs(diff(signal))>200];
    badIntervals = FindInterval(bad); badIntervals = [badIntervals(:,1)-1 badIntervals(:,2)+1];
    badIntervals = ConsolidateIntervalsFast(badIntervals,'epsilon',15);
    save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
    datestr((datenum(clock)))
    noiseIntervalIndices = badIntervals;
    noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;
    m = memmapfile(datFile, 'Format','int16','Writable',true);
    for j = 1:nChannels
        badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
        goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
        badIndices = sub2ind([nChannels,nSamples],j*ones(size(badTimeIndices)),badTimeIndices);
        goodIndices = sub2ind([nChannels,nSamples],j*ones(size(goodTimeIndices)),goodTimeIndices);
        goodValues = m.Data(goodIndices);
        interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
        m.Data(badIndices) = int16(interpolated);
    end
    cleaned(i,2) = true;
    disp(['dat-file for ' basename ' cleaned! Starting a new Kilosort']);
    catch
        disp([datestr(clock) ': Problem with cleaning dat-file for ' basename]);
    end
end

%% Kilosort:

for i=41:length(folders)
    try
    basepath = folders{i};
    cd(basepath)
    f = dir('Kilosort*');
    if ~isempty(f),done(i,1) = true; kilosorted(i,1) = true; continue; end
    cd(basepath);
    basepath = folders{i};
    [~,basename] = fileparts(basepath);
    disp([datestr(clock) ': Starting kilosort for ' basename]);
    session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
    save(fullfile(basepath,[basename '.session.mat']),'session')
    % Kilosort concatenated sessions
    if isempty(dir('KilosortGT*')) && spikeSort
        kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
        load(fullfile(kilosortFolder,'rez.mat'),'rez');
        CleanRez(rez,'savepath',kilosortFolder);
        %     PhyAutoClustering(kilosortFolder);
    end
    kilosorted(i,1) = true;
    disp([datestr(clock) ': Kilosort done for ' basename]);

    % Get brain states
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
    done(i,1) = true;
    catch
        disp('Problem with session...');
    end
end

%% Old code
% nDays = 54;
% 
% folders = cell(nDays,1);
% for i=1:length(folders),
%     folders{i} = ['N:\V1test\AO52\day' num2str(i)];
% %     folders{i} = ['N:\V1test\V1Jean\day' num2str(i)];
% end
% 
% fillMissingDatFiles = true;
% fillTypes = [];
% analogInputs = true;
% analogChannels = 1;
% digitalInputs = true;
% digitalChannels = 1:4;
% getAcceleration = true;
% cleanArtifacts = false;
% stateScore = true;
% spikeSort = true;
% getPos = false;
% removeNoise = false;
% runSummary = false;
% SSD_path = 'D:\KiloSort';
% done = false(nDays,1);
% for i=1:nDays
%     try
%         basepath = folders{i}
%         basename = basenameFromBasepath(basepath);
% 
%         cd(basepath);
%         f = dir('Kilosort*');
%         if isempty(f) % No kilosort done, launch all preprocessing:
%             disp(['Launching preprossing for ' basename]);
%             %% Pull meta data
%             % Get session names
%             if strcmp(basepath(end),filesep)
%                 basepath = basepath(1:end-1);
%             end
%             [~,basename] = fileparts(basepath);
% 
%             % Get xml file from the parent folder
%             [parentFolder,dayName] = fileparts(basepath);
%             [~,projectName] = fileparts(parentFolder);
%             parentXml = fullfile(parentFolder,[projectName '.xml']);
%             xmlFile = fullfile(basepath,[dayName '.xml']);
%             if ~exist(xmlFile,'file'), copyfile(parentXml,xmlFile); end
%             rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
%             rhdFile = rhdFile(1);
%             if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
%                 copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
%             end
% 
%             session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
%             save([basename '.session.mat'],'session');
%             %% Fill missing dat files of zeros
%             if fillMissingDatFiles
%                 if isempty(fillTypes)
%                     fillTypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
%                 end
%                 for ii = 1:length(fillTypes)
%                     fillMissingDats('basepath',basepath,'fileType',fillTypes{ii});
%                 end
%             end
%             %% Concatenate sessions
%             cd(basepath);
% 
%             disp('Concatenate session folders...');
%             concatenateDats(pwd,0,0);
% 
%             % Process additional inputs - CHECK FOR OUR LAB
% 
%             if analogInputs
%                 try
%                     if  ~isempty(analogChannels)
%                         analogInp = computeAnalogInputs('analogCh',analogChannels,'saveMat',true,'fs',session.extracellular.sr);
%                     else
%                         analogInp = computeAnalogInputs('analogCh',[],'saveMat',true,'fs',session.extracellular.sr);
%                     end
% 
%                     % analog pulses ...
%                     [pulses] = getAnalogPulses('samplingRate',session.extracellular.sr);
%                 end
%             end
% 
%             % Digital inputs
%             try
%                 digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
%             catch
%                 display('Error with digital inputs. This step was skipped.');
%             end
%             % Auxilary input
%             if getAcceleration
%                 try
%                     accel = computeIntanAccel('saveMat',true);
%                 end
%             end
% 
%             % Make LFP
%             LFPfromDat(pwd,'outFs',1250,'useGPU',false);
% 
%             %% Clean data  - CHECK FOR OUR LAB
%             % Remove stimulation artifacts
%             if cleanArtifacts && analogInputs
%                 [pulses] = getAnalogPulses(analogInp,'analogCh',analogChannels);
%                 cleanPulses(pulses.ints{1}(:));
%             end
% 
%             % remove noise from data for cleaner spike sorting
%             % if removeNoise
%             %     NoiseRemoval(pwd); % not very well tested yet: Raly: this is bad
%             % end
% 
%         else % Kilosort has been done before, check if session has already been cleaned:
%             predone(i,1) = true;
%             f_noise = dir('noiseIntervalsDat.events.mat');
%             if ~isempty(f_noise)
%                 disp([basename ' has already been cleaned']);
%                 done(i,1) = true;
%                 error('done'); % continue to next session (we are in a try/catch statement)
%             else
%                 disp(['Renaming the old kilosort of ' basename]);
%                 % Move all the previous results for posterity:
%                 movefile(f.name,strrep(f.name,'Kilosort','raw_Kilosort'));
%                 % if there are any CellExplorer files, move them to that folder:
%                 x = dir('*session.mat');
%                 y = dir('*.cellinfo.mat');
%                 files2move = [x;y];
%                 if ~isempty(files2move)
%                     disp(['Moving cellinfo files of ' basename]);
%                     rawpath = fullfile(basepath,strrep(f.name,'Kilosort','raw_Kilosort'));
%                     cd(rawpath);
%                     mkdir('CellExplorer');
%                     for k=1:size(files2move,1)
%                         movefile(fullfile(basepath,files2move(k).name),fullfile(rawpath,'CellExplorer',files2move(k).name));
%                     end
%                 end
%             end
%         end
%         disp(['Cleaning dat-file for ' basename]);
%         % Apply AO52-specific clean:
%         rejectChannels = 1+[8 12 14 15]; % for AO52
%         nChannels = 64;
% 
%         basename = basenameFromBasepath(basepath);
%         datFile = [basepath,filesep, basename, '.dat'];
%         m = memmapfile(datFile, 'Format','int16','Writable',true);
%         data = reshape(m.data,64,[]);
%         nSamples = size(data,2);
%         okChannels = ~ismember((1:size(data,1))',rejectChannels);
%         signal = mean(data(okChannels,:))';
%         bad = [false; abs(diff(signal))>200];
%         badIntervals = FindInterval(bad); badIntervals = [badIntervals(:,1)-1 badIntervals(:,2)+1];
%         badIntervals = ConsolidateIntervalsFast(badIntervals,'epsilon',15);
%         save(fullfile(basepath,'noiseIntervalsDat.events.mat'),'badIntervals');
%         datestr((datenum(clock)))
%         noiseIntervalIndices = badIntervals;
%         noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;
%         m = memmapfile(datFile, 'Format','int16','Writable',true);
%         for j = 1:nChannels
%             badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
%             goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
%             badIndices = sub2ind([nChannels,nSamples],j*ones(size(badTimeIndices)),badTimeIndices);
%             goodIndices = sub2ind([nChannels,nSamples],j*ones(size(goodTimeIndices)),goodTimeIndices);
%             goodValues = m.Data(goodIndices);
%             interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
%             m.Data(badIndices) = int16(interpolated);
%         end
%         cleaned(i,1) = true;
%         disp(['dat-file for ' basename ' cleaned! Starting a new Kilosort']);
%         
%         session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
%         save(fullfile(basepath,[basenameFromBasepath(basepath) '.session.mat']),'session')
%         % Kilosort concatenated sessions
%         if isempty(dir('KilosortGT*')) && spikeSort
%             kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
%             load(fullfile(kilosortFolder,'rez.mat'),'rez');
%             CleanRez(rez,'savepath',kilosortFolder);
%             %     PhyAutoClustering(kilosortFolder);
%         end
%         kilosorted(i,1) = true;
% 
%         % Get brain states
%         % an automatic way of flaging bad channels is needed
%         if stateScore
%             try
%                 if exist('pulses','var')
%                     SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
%                     thetaEpochs(pwd);
%                 else
%                     SleepScoreMaster(pwd,'noPrompts',true); % takes lfp in base 0
%                     thetaEpochs(pwd);
%                 end
%             catch
%                 disp('Problem with SleepScore skyping...');
%             end
%         end
% 
%         done(i,1) = true;
%     catch
%         warning(['Issue with session ' num2str(i)]);
%     end
% end
