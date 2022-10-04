function reconcatenateDats(basepath,sortFiles,otherdattypes,SSD_path)
% reconcatenateDats - Function to help you concatenate auxiliary .dat files that were missed during preprocessing
%
% Sometimes after preprocessing a session, you may realize that there was 
% an input that was not concatenated. To avoid running concatenateDats from 
% scratch (as it takes a long time), run concatenateOtherDats to fill in
% the missing input .dat files. 
%
%  USAGE
%
%    concatenateDats(basepath,sortFiles)
%
%  INPUTS
%
%    basepath                  path to session folder
%    sortFiles                 boolean denoting whether to sort files according 
%                              to time of recording (1) or not (0) and thus sort 
%                              them alphabetically.
%    otherdattypes             a cell containing the filename that this function
%                              should look for in other subfolders. See examples below.
%                              
%                              
%
%  OUTPUT
%     Concatenated files will be saved in the basepath folder.
%
%  EXAMPLES
%
%     % To concatenate digitalin.dat and analogin.dat files from individual subfolders: (e.g. after running "fillMissingDats")
%     concatenateOtherDats(basepath,false,{'digitalin','analogin'}); 
%
%     % To (re)concatenate the recording .dat (e.g. if it has been corrupted/modified/deleted)
%     concatenateOtherDats(basepath,false,{'amplifier'}); 
%     % Don't forget to rename it once it's concatenated as a "amplifier.dat" in the basepath folder
%
%
% Copyright (C) 2022 Ralitsa Todorova (modified from concatenateDats)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('SSD_path','var')
    SSD = true; else, SSD = false;
end

bad_otherdattypes = [];
if SSD
    for odidx = 1:length(otherdattypes)
        %eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
        newpaths.(otherdattypes{odidx}) = fullfile(SSD_path,[basenameFromBasepath(basepath) '_' otherdattypes{odidx}, '.dat']);
        newpathsFinal.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
    end
else
    for odidx = 1:length(otherdattypes)
        %eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
        newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
    end
end


d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders

for a = 1:length(d)
    %look in each subfolder
    if any(~ismember(d(a).name,'.')) && d(a).isdir
        %Check for amplifier.dat or subfolderbaseName.dat
        if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
            ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
        elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
            ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
        else
            ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
            digitalinfile = fullfile(basepath,d(a).name,'digitalin.dat');
        end
        if ~exist(ampfile,'file') && exist(digitalinfile,'file')
            % create an empty ampfile. This will allow us to preprocess sessions without recorded brain activity (but e.g. digital events)
            if ~exist(ampfile,'file'), fclose(fopen(ampfile, 'w')); end
        end
        if exist(ampfile,'file')
            rcount = rcount+1;
            datpaths.amplifier{rcount} = ampfile;
            t = dir(ampfile);
            datsizes.amplifier(rcount) = t.bytes;
            recordingnames{rcount} = d(a).name;
            for odidx = 1:length(otherdattypes)%loop through other .dat types found here
                % eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
                datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
                %eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
                d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
                if isempty(d2)
                    bad_otherdattypes(odidx) = 1;
                else
                    %eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
                    datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
                end
            end
        end
    end
end

%%
if sortFiles
    try
        names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
        names2sort = cell2mat(names2sort);
%         if isempty(names2sort{1}) && ~isempty(recordingnames{1})
%             error('Last 6 digits were not numeric and therefore do not reflect the recording time.');
%         end
        disp('Assuming the last 6 digits reflect recording time.')
        %disp('Don''t like it? Write in some new options for sorting.')
    catch
       % names2sort = 1:length(recordingnames);
        disp('Last 6 digits not numeric... sorting alphabetically')
    end

    [~,I] = sort(names2sort);
    recordingnames = recordingnames(I);
    for odidx = 1:length(otherdattypes)
        datpaths.(otherdattypes{odidx}) = datpaths.(otherdattypes{odidx})(I);
        datsizes.(otherdattypes{odidx}) = datsizes.(otherdattypes{odidx})(I);
    end
end

disp('Concatenating Other Dats..... continue to be patient')
for odidx = 1:length(otherdattypes)
    if isunix
        cs = strjoin(datpaths.(otherdattypes{odidx}));
        catstring = ['! cat ', cs, ' > ',newpaths.(otherdattypes{odidx})];
    elseif ispc
        if length(datpaths.(otherdattypes{odidx}))>1
            for didx = 1:length(datpaths.(otherdattypes{odidx}))
                datpathsplus{didx} = [datpaths.(otherdattypes{odidx}){didx} '+'];
            end
            %Last file string shouldn't end with '+'
            datpathsplus{length(datpaths.(otherdattypes{odidx}))} = datpaths.(otherdattypes{odidx}){didx};
        else
            datpathsplus = datpaths.(otherdattypes{odidx});
        end
        cs = strjoin(datpathsplus);
        catstring = ['! copy /b ', cs, ' ',newpaths.(otherdattypes{odidx})];
    end
    eval(catstring)%execute concatenation
    % Check that size of resultant .dat is equal to the sum of the components
    t = dir(newpaths.(otherdattypes{odidx}));
    if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
        error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
        deleteoriginaldatsbool = 0;
        sizecheck.(otherdattypes{odidx}) = false;
    else
        sizecheck.(otherdattypes{odidx}) = true;
        disp([otherdattypes{odidx} ' concatenated and size checked'])
    end
end

if SSD
    for odidx = 1:length(otherdattypes)
        disp([datestr(clock) '. Copying local dat file to ' basepath '...']);
        copyfile(newpaths.(otherdattypes{odidx}),newpathsFinal.(otherdattypes{odidx}));
    end
end

