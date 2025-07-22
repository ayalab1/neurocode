% first run the following to merge rhd files:
% neurocode\preProcessing\intan\Windows\IntanFileMerger.exe
% run read_Intan_RHD2000_file.m to load data into memory
% save all data file using the following

% files_table = table();
% files_table.files = {'amplifier', 'auxiliary', 'digitalin', 'digitalout', 'analogin', 'time', 'supply'}';
% files_table.data_type = {'int16', 'uint16', 'uint16', 'uint16', 'uint16', 'int32', 'uint16'}';

f = fopen('amplifier.dat','w'); 
fwrite(f,amplifier_data,'int16');
fclose(f);

f = fopen('auxiliary.dat','w'); 
fwrite(f,aux_input_data(:),'uint16');
fclose(f);

f = fopen('analogin.dat','w'); 
fwrite(f,board_adc_data(:),'uint16');
fclose(f);

f = fopen('digitalin.dat','w'); 
fwrite(f,board_dig_in_data(:),'uint16');
fclose(f);

f = fopen('time.dat','w'); 
fwrite(f,t_amplifier(:),'int32');
fclose(f);
