function cutDatChans(fnameIn,fnameOut,nChans,chans)
% cutDatChans(fnameIn,fnameOut,nChans,chans)
% create a new dat file with a subset of channels (chans). 
% Call it fromt the folder where the dat file is.

% THIS SHOULD BE CHANGED TO TO READ ADN WRITE DAT FILE IN CHUNKS

file = memmapfile(fnameIn,'Format','int16','Writable',false);
datain = reshape(file.data,nChans,[]);

dataout = datain(chans,:);

f = fopen(fnameOut,'w'); 
fwrite(f,dataout(:),'int16');
fclose(f);


end