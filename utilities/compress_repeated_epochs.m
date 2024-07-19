function results = compress_repeated_epochs(epoch_df, epoch_name)
% compress_repeated_epochs: Compresses epoch_df
% If there are back to back epochs of the same name, it will combine them
%
% Input:
%   epoch_df (table): MATLAB table with epoch data, containing at least the following columns:
%       - environment (cell array of strings): The name of the environment for each epoch
%       - startTime (numeric array): The start time of each epoch
%       - stopTime (numeric array): The stop time of each epoch
%   epoch_name (optional, string): Specific name of the epoch environment to compress. If not provided, all environments are considered.
%
% Output:
%   results (table): Compressed epoch data with combined back-to-back epochs
%
% Example:
%   epoch_df = load_epoch()
%                      name                      startTime    stopTime      environment      behavioralParadigm     manipulation        stimuli               notes                                basepath
% _______________________________________    _________    ________    _______________    __________________    _______________    ____________    _____________________    ________________________________________________
%
% {'test_DC_pulse_1V_20mA_240506_105224'}          0       546.23     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'test_DC_pulse_1V_15mA_240506_110149'}     546.23       654.89     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'test_DC_pulse_1V_17mA_240506_110431'}     654.89        794.2     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'test_DC_pulse_1V_19mA_240506_110731'}      794.2       954.44     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'presleep_240506_111315'             }     954.44       9096.1     {'sleep'      }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'cheeseboard_240506_133325'          }     9096.1        10444     {'cheeseboard'}       {'25'      }       {'pfc_silence'}    {0×0 double}    {'square_wave_pulse'}    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'cheeseboard_2_240506_135824'        }      10444        11457     {'cheeseboard'}       {'25'      }       {'pfc_silence'}    {0×0 double}    {'square_wave_pulse'}    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'postsleep_240506_141807'            }      11457        20683     {'sleep'      }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
%
%   results = compress_repeated_epochs(epoch_df)
%
%                  name                      startTime    stopTime      environment      behavioralParadigm     manipulation        stimuli               notes                                basepath
% _______________________________________    _________    ________    _______________    __________________    _______________    ____________    _____________________    ________________________________________________
%
% {'test_DC_pulse_1V_20mA_240506_105224'}          0       954.44     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'presleep_240506_111315'             }     954.44       9096.1     {'sleep'      }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'cheeseboard_240506_133325'          }     9096.1        11457     {'cheeseboard'}       {'25'      }       {'pfc_silence'}    {0×0 double}    {'square_wave_pulse'}    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'postsleep_240506_141807'            }      11457        20683     {'sleep'      }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
%
%
%
%   results = compress_repeated_epochs(epoch_df, 'sleep')
%                      name                      startTime    stopTime      environment      behavioralParadigm     manipulation        stimuli               notes                                basepath
% _______________________________________    _________    ________    _______________    __________________    _______________    ____________    _____________________    ________________________________________________
%
% {'test_DC_pulse_1V_20mA_240506_105224'}          0       546.23     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'test_DC_pulse_1V_15mA_240506_110149'}     546.23       654.89     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'test_DC_pulse_1V_17mA_240506_110431'}     654.89        794.2     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'test_DC_pulse_1V_19mA_240506_110731'}      794.2       954.44     {'homecage'   }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'presleep_240506_111315'             }     954.44       9096.1     {'sleep'      }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'cheeseboard_240506_133325'          }     9096.1        10444     {'cheeseboard'}       {'25'      }       {'pfc_silence'}    {0×0 double}    {'square_wave_pulse'}    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'cheeseboard_2_240506_135824'        }      10444        11457     {'cheeseboard'}       {'25'      }       {'pfc_silence'}    {0×0 double}    {'square_wave_pulse'}    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
% {'postsleep_240506_141807'            }      11457        20683     {'sleep'      }       {0×0 double}       {0×0 cell     }    {0×0 double}    {0×0 double         }    {'U:\data\hpc_ctx_project\HP04\day_36_20240506'}
%
% adapted from neuro_py.session.locate_epoch.compress_repeated_epochs.py
% Ryan H, Laura B

% Check if epoch_name is provided, otherwise set it to empty
if nargin < 2
    epoch_name = [];
end

% Initialize match array with NaNs
match = NaN(size(epoch_df.environment, 1), 1);

% Loop through the epochs to find consecutive matches
for i = 1:height(epoch_df) - 1
    if isnan(match(i))
        if isempty(epoch_name) % If no specific epoch_name is provided
            if strcmp(epoch_df.environment{i}, epoch_df.environment{i+1})
                match(i) = i;
                match(i+1) = i;
                for match_i = i + 2:height(epoch_df)
                    if match_i > height(epoch_df) || ~strcmp(epoch_df.environment{i}, epoch_df.environment{match_i})
                        break;
                    end
                    match(match_i) = i;
                end
            end
        else % If a specific epoch_name is provided
            if strcmp(epoch_df.environment{i}, epoch_df.environment{i+1}) && strcmp(epoch_df.environment{i}, epoch_name)
                match(i) = i;
                match(i+1) = i;
                for match_i = i + 2:height(epoch_df)
                    if match_i > height(epoch_df) || ~strcmp(epoch_df.environment{i}, epoch_df.environment{match_i})
                        break;
                    end
                    match(match_i) = i;
                end
            end
        end
    end
end

% Replace remaining NaNs in match array with unique large numbers
match(isnan(match)) = (1:length(match(isnan(match))))' * 2000;

% Initialize results table
results = table();
no_nan_match = match(~isnan(match));
unique_matches = unique(no_nan_match, 'stable'); % Ensure stable sorting

% Combine epochs based on matches
for m = 1:length(unique_matches)
    temp_dict = struct();
    idx = match == unique_matches(m);
    vars = epoch_df.Properties.VariableNames;

    % Collect data for the combined epoch
    for item = 1:length(vars)
        temp_dict.(vars{item}) = epoch_df.(vars{item})(find(idx, 1));
    end

    % Set startTime to the earliest startTime and stopTime to the latest stopTime
    temp_dict.startTime = min(epoch_df.startTime(idx));
    temp_dict.stopTime = max(epoch_df.stopTime(idx));

    % Convert the struct to a table
    temp_table = struct2table(temp_dict);

    % Append the combined epoch to the results table
    if isempty(results)
        results = temp_table;
    else
        results = [results; temp_table];
    end
end
end
