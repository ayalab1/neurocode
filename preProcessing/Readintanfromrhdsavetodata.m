% step1 run
% D:\Documents\GitHub\neurocode\preProcessing\intan\Windows\IntanFileMerger.exe
% to merge all rhd 
% step1 run read_Intan_RHD2000_file.m
%save all data file
f = fopen('amplifier.dat','w'); % new file should be test.dat
fwrite(f,amplifier_data,'int16');
fclose(f);
f = fopen('auxiliary.dat','w'); % new file should be test.dat
fwrite(f,aux_input_data,'int16');
fclose(f);





%% for digital in directly save as digital events file
% board_dig_in_data(11:16,1:length(board_dig_in_data))=0;
% f = fopen('digitalin2.dat','w'); % new file should be test.dat
% fwrite(f,board_dig_in_data,'int16');
% fclose(f);
digitalChannels=[1,2,3,4,5,6,7,8,9,10];
digitalInp = getDigitalIn_hy(digitalChannels,'fs',20000);

%individually detect digital events file
mainFolder = 'X:/data/PPP/PPP9/day7';
items = dir(mainFolder) ;

    % Iterate through each item in the folder
    for i = 1:length(items)
        item = items(i);
        
        % Ignore the current and parent folder references
        if strcmp(item.name, '.') || strcmp(item.name, '..')
            continue;
        end       
        % Create the full path of the current item
        itemPath = fullfile(Folder , item.name);
        
        % Check if the current item is a folder
        if item.isdir
            % Recursively call the folder processing code to process the subfolder

        % Check if the current item is a .dat file
        elseif strcmpi(item.name(end-3:end), '.dat')
            % Run your desired function on the .dat file
            % Replace 'YourFunction' with the actual function name and pass the itemPath as an argument
            cd(itemPath)
            digitalChannels=[1,2,3,4,5,6,7,8,9,10,11,12];
            digitalInp = getDigitalIn_hy(digitalChannels,'fs',20000);;
        end
    end

%merge all digital events file 
basepath=pwd;
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.MergePoints.events.mat']))
sumT=0;
digitalInp= cell(1, 11);
digitalduration= cell(1, 11);

for ii=1:size(MergePoints.foldernames ,2)
    load([basepath filesep MergePoints.foldernames{1,ii} filesep 'digitalIn.events.mat'])
    %oprerate each cell
    sumT=MergePoints.timestamps(ii,1);
    for cellIndex = 1:numel(digitalIn.timestampsOn)
        digitalIn.timestampsOn{1,cellIndex} = digitalIn.timestampsOn{1,cellIndex} + sumT;
        digitalInp{1, cellIndex}=[digitalInp{1, cellIndex};digitalIn.timestampsOn{1,cellIndex}];
        digitalduration{1, cellIndex}=[digitalduration{1, cellIndex} digitalIn.dur{1,cellIndex}];
    end
    clear digitalIn
end
digitalIn.timestampsOn=digitalInp;
digitalIn.duration=digitalduration;
save('digitalIn.events.mat','digitalIn')



% dealling with time.dat
%timeduration from amplifier file 
% step1: Get num_sample from concatenateDats
%timestamps_samples=numsamples-1
%duration=timestamps_samples/20000
% step2:
MergePoints.timestamps(2,2)=MergePoints.timestamps(2,1)+ -5.0000e-05+numsamples(1,2)/20000
for i=3:7
    for ii=1:2
        MergePoints.timestamps(i,ii)=MergePoints.timestamps(i,ii)+ -5.0000e-05+numsamples(1,2)/20000
    end
end

MergePoints.timestamps_samples(2,2)=MergePoints.timestamps_samples(2,1)+numsamples(1,2)-1
for i=3:7
    for ii=1:2
        MergePoints.timestamps_samples(i,ii)=MergePoints.timestamps_samples(i,ii)+numsamples(1,2)-1
    end
end

f = fopen('time.dat','w'); % new file should be test.dat
fwrite(f,t_amplifier,'int32');
fclose(f);



f = fopen('analogin.dat','w'); % new file should be test.dat
fwrite(f,data(:),'int16');
fclose(f);