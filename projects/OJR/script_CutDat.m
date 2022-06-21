
nChannels = ?; % set this manually
maxTime = ?; % set manually: this should be in seconds

% rename files: do this only once as doing it more than once would remove the original files!!!
if exist('amplifier_original.dat','file')
    error('amplifier_original.dat file already exists. Don''t run script again, it might overwrite your original files!!!')
else
    movefile('auxiliary.dat','auxiliary_original.dat');
    movefile('digitalin.dat','digitalin_original.dat'); 
    movefile('analogin.dat','analogin_original.dat');
    movefile('time.dat','time_original.dat');
    movefile('amplifier.dat','amplifier_original.dat');
end


file = memmapfile('amplifier_original.dat','Format','int16','Writable',false);
data = reshape(file.data,nChannels,[]);
originalSize = size(data,2);
data = data(:,1:ceil(maxTime*20000));
f = fopen('amplifier.dat','w'); % new file should be test.dat
fwrite(f,data(:),'int16');
fclose(f);

file = memmapfile('auxiliary_original.dat','Format','int16','Writable',false);
nChannels = length(file.data)/originalSize;
data = reshape(file.data,nChannels,[]);
data = data(:,1:ceil(maxTime*20000));
f = fopen('auxiliary.dat','w'); % new file should be test.dat
fwrite(f,data(:),'int16');
fclose(f);

file = memmapfile('digitalin_original.dat','Format','int16','Writable',false);
nChannels = length(file.data)/originalSize;
data = reshape(file.data,nChannels,[]);
data = data(:,1:ceil(maxTime*20000));
f = fopen('digitalin.dat','w'); % new file should be test.dat
fwrite(f,data(:),'int16');
fclose(f);

file = memmapfile('time_original.dat','Format','int16','Writable',false);
nChannels = length(file.data)/originalSize;
data = reshape(file.data,nChannels,[]);
data = data(:,1:ceil(maxTime*20000));
f = fopen('time.dat','w'); % new file should be test.dat
fwrite(f,data(:),'int16');
fclose(f);

file = memmapfile('analogin_original.dat','Format','int16','Writable',false);
nChannels = length(file.data)/originalSize;
data = reshape(file.data,nChannels,[]);
data = data(:,1:ceil(maxTime*20000));
f = fopen('analogin.dat','w'); % new file should be test.dat
fwrite(f,data(:),'int16');
fclose(f);
