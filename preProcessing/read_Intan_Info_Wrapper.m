function varargout = read_Intan_Info_Wrapper(filename, varargin)

p = inputParser;

addParameter(p,'returnStruct',true,@islogical)
parse(p,varargin{:});
returnStruct = p.Results.returnStruct;

try % temporary solution, need to figure out whats the problem
[amplifier_channels, notes, aux_input_channels, spike_triggers,...
board_dig_in_channels, board_adc_channels, supply_voltage_channels, frequency_parameters ] =...    
read_Intan_RHD2000_file([],filename);
catch 
[amplifier_channels, notes, aux_input_channels, spike_triggers,...
board_dig_in_channels, board_adc_channels, supply_voltage_channels, frequency_parameters ] =...    
read_Intan_RHD2000_file_bz;
end

if returnStruct
    data.amplifier_channels = amplifier_channels;
    data.notes = notes;
    data.aux_input_channels = aux_input_channels;
    data.spike_triggers = spike_triggers;
    data.board_dig_in_channels = board_dig_in_channels;
    data.board_adc_channels = board_adc_channels;
    data.supply_voltage_channels = supply_voltage_channels;
    data.frequency_parameters = frequency_parameters;
    
    varargout = cell(1,1);
    varargout{1} = data;
else
    varargout = cell(1,8);
    varargout{1} = amplifier_channels;
    varargout{2} = notes;
    varargout{3} = aux_input_channels;
    varargout{4} = spike_triggers;
    varargout{5} = board_dig_in_channels;
    varargout{6} = board_adc_channels;
    varargout{7} = supply_voltage_channels;
    varargout{8} = frequency_parameters;
end

    
    