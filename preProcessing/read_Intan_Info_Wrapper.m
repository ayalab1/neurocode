function varargout = read_Intan_Info_Wrapper(filename, varargin)
%
%
% [varargout] = [read_Intan_Info_Wrapper(filename, varargin)]
%  
%
% [This is a wrapper for read_Intan_RHD2000_file_bz output. The output
% should be a structure that can be read out by this function]
%
%  USAGE
%    varargout = cell(1,8);
%    [amplifier_channels, notes, aux_input_channels, spike_triggers,...
%    board_dig_in_channels, supply_voltage_channels, frequency_parameters,board_adc_channels ]
%     = read_Intan_Info_Wrapper(returnstruct)]
%
%
%  INPUTS
%    [parser]      [parser, see options list]
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     ['returnstruct']  [this is a struct from read_Intan_RHD2000_file_bz
%                        output that this function reads]
%    =========================================================================
% 
%  OUTPUTS
%    [amplifier_channels]     [amplifier channels]
%    [notes]                  [notes]
%    [aux_input_channels]     [auxiliary input channels]
%    [spike_triggers]         [spike triggers]
%    [board_dig_in_channels]  [DIG's from board as channels]
%    [supply_voltage_channels][voltage supply channels]
%    [frequency_parameters]   [output of frequency parameters]
%    [board_adc_channels]     [output of INTAN adc as channels]
%
%  NOTE
%  [this is a wrapper to use read_Intan_RHD2000_file_bz]
%
%  SEE ALSO
%
% [AntonioFR] [2021- 2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%-------------------------------------------------------------------------


p = inputParser;

addParameter(p,'returnStruct',true,@islogical)
parse(p,varargin{:});
returnStruct = p.Results.returnStruct;

[amplifier_channels, notes, aux_input_channels, spike_triggers,...
board_dig_in_channels, supply_voltage_channels, frequency_parameters,board_adc_channels] =...    
read_Intan_RHD2000_file_bz;


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

    
    