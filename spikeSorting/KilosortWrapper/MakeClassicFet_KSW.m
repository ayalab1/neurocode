function MakeClassicFet_KSW(basename,dirname)

%get basic input params if not there
if ~exist('basename','var')
    [~,basename] = fileparts(cd);
end
if ~exist('dirname','var')
    dirname = cd;
end

cd(dirname)

% if already .fets, move them to a new folder called "PreviousFets"
d = dir(fullfile(dirname,[basename '.fet.*']));
if ~isempty(d);
    mkdir(fullfile(dirname,'PreviousFets'))
    for a = 1:length(d);
        movefile(fullfile(dirname,d(a).name),fullfile(dirname,'PreviousFets',d(a).name));
    end
end

%if ~exist([basename '.fil'],'file')
%    cmd = ['process_mhipass ' basename];
%    system(cmd)
%end

if ~exist([basename '.fil'])
    xml =  LoadXml([basename '.xml']);
    inname = [basename '.dat'];
    outname = [basename '.fil'];
    numch = num2str(xml.nChannels);
    sampl = num2str(xml.SampleRate);
    lowband = '800';
    highband = '9500';
    forder = '50';
    gain = '1';
    offset = '0';
    
    firfilter_KSW(inname,outname,numch,sampl,lowband,highband,forder,gain,offset);
end    


%get xml
[xml, rxml] = LoadXml([basename '.xml']);
nGrps = length(xml.SpkGrps);
for g = 1:nGrps
    cmd = ['process_pca_multi ' basename ' ' num2str(g)];
%     cmd = ['ndm_pca ' basename ' ' num2str(i)];
    disp(cmd)
    a =  system(cmd)
end


if a ==0
    cmd = ['rm ' basename '.fil'];
    system(cmd)
end

