function concatenateOtherDats(basepath);

deleteoriginaldatsbool = 0;
sortFiles = 0;
basename = basenameFromBasepath(basepath);
otherdattypes = {'analogin';'digitalin';'auxiliary';'time'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
    %eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
    newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
newpaths
bad_otherdattypes
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

disp('Concatenating Other Dats..... continue to be patient')
for odidx = 1:length(otherdattypes)
    if isunix
        cs = strjoin(datpaths.(otherdattypes{odidx}));
        catstring = ['! cat ', cs, ' > ',newpaths.(otherdattypes{odidx})];
    elseif ispc%As of 4/9/2017 - never tested
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
