function fname = KiloSortLinuxDir(basename,basepath,gpuDeviceNum)


%get size of dat
[a,b] = system(['stat --format=%s "' basepath '/' basename '.dat"']);
datsize = str2num(b)/1e9;

%get all hard drives

tmp = dir('/sys/block');
kp = arrayfun(@(a) any(regexp(a.name,'sd')),tmp);
HDs = tmp(kp);

%find all SSDs
for i = 1:length(HDs)
    [a,b] = system(['cat /sys/block/' HDs(i).name '/queue/rotational']);
    isSSD(i) = ~(str2num(b));
    
end





%get free space on all mounts
freespace = [];
idx = 1;
for i = 1:length(HDs)
    
    HD = dir('/dev/');
    kp = arrayfun(@(a) any(regexp(a.name,[HDs(i).name '[0-9]'])),HD);
    HD = HD(kp);
    
    
    
    for j = 1:length(HD)
        [a,mount] =  system(['grep "/dev/' HD(j).name '" /proc/mounts | cut -d '' '' -f 2']);
        if ~any(regexpi(mount,'/.')) && ~isempty(mount)
            %store on non root SSD
            
            [a,b] = system(['df -Ph ' mount ' | tail -1 | awk ''{print $4}''']);
            freespace(idx) = str2num(b(regexp(b,'[0-9]')));
            mnt{idx} = mount;
            SSD(idx) = isSSD(i);
            idx = idx+1;
        elseif any(regexpi(mount,'/.'))
            [a,b] = system('df -Ph /home | tail -1 | awk ''{print $4}''');
            freespace(idx) = str2num(b(regexp(b,'[0-9]')));
            SSD(idx) = isSSD(i);
            [a,b] = system('whoami');
            username = b(1:regexp(b,'\n')-1);
            
            
            mnt{idx} = ['/home/' username];
            idx = idx+1;
        end
        
    end
end

%priorize SSDs

mountSSD = mnt(SSD);
freespaceSSD = freespace(SSD);

mountHD = mnt(~SSD);
freespaceHD = freespace(~SSD);


if any( (freespaceSSD-datsize) > .5)
    %save 500MB on the SSD, can be decreased
    [~,b] = max(freespaceSSD-datsize);
    
    fname = [mountSSD{b} '/temp_wh_' num2str(gpuDeviceNum) '.dat'];
elseif  any( (freespaceHD-datsize) > .5)
    [~,b] = max(freespaceSSD-datsize);
    
    fname = [mountHD{b} '/temp_wh_' num2str(gpuDeviceNum) '.dat'];
else
    error('NO DISK SPACE')
    
end






end