function batch_analysis(df, save_path, func, varargin)
% batch_analysis: run the main loop in parallel or serial
% Inputs:
%    df: table with basepath column
%    save_path: str path to save results
%    func: function to run on each basepath in df
%    varargin: additional arguments to configure the function and to pass to func
%
% example:
%
% Select sessions to run. sessions is a table with column "basepath"
% sessions = readtable("Z:\home\ryanh\projects\hpc_ctx\opto_sessions.csv", "Delimiter", ',')
%
% Select path to save results. Project specific path so to not pollute data paths
% save_path = 'Z:\home\ryanh\projects\hpc_ctx\theta_seq\theta_seq_v1';
%
% run analysis on all sessions. Here the analysis is "run" and takes basepath and other varargin
% batch_analysis(sessions, save_path, @run)

% Create input parser
p = inputParser();
p.KeepUnmatched = true;
% Add required arguments
addRequired(p, 'df', @(x) istable(x) && any(strcmp('basepath', x.Properties.VariableNames)));
addRequired(p, 'save_path', @ischar);
addRequired(p, 'func', @(x) isa(x, 'function_handle'));

% Add optional arguments with default values
addParameter(p, 'parallel', true, @islogical);
addParameter(p, 'verbose', false, @islogical);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'skip_if_error', false, @islogical);
addParameter(p, 'num_cores', [], @(x) isempty(x) || isnumeric(x));

% Parse inputs
parse(p, df, save_path, func, varargin{:});

% Extract values
df = p.Results.df;
save_path = p.Results.save_path;
func = p.Results.func;
parallel = p.Results.parallel;
verbose = p.Results.verbose;
overwrite = p.Results.overwrite;
skip_if_error = p.Results.skip_if_error;
num_cores = p.Results.num_cores;

% Remaining arguments are for func
% struct to cell
varargin = {};
vars = fields(p.Unmatched);
for field_i = 1:length(vars)
    varargin = [varargin, vars{field_i}];
    varargin = [varargin, p.Unmatched.(vars{field_i})];
end

% Find sessions to run
basepaths = unique(df.basepath);
% Create save_path if it doesn't exist
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Run in parallel if parallel is True
if parallel
    % Get number of cores
    if isempty(num_cores)
        num_cores = feature('numcores');
    end
    WaitMessage = parfor_wait(length(basepaths));
    % Run in parallel
    parfor (i = 1:length(basepaths), num_cores)
        main_loop(basepaths{i}, save_path, func, overwrite, skip_if_error, varargin{:});
        WaitMessage.Send;
    end
else
    if verbose
        WaitMessage = parfor_wait(length(basepaths));
    end
    % Run in serial
    for i = 1:length(basepaths)
        if verbose
            disp(basepaths{i});
        end
        % Run main_loop on each basepath in df
        main_loop(basepaths{i}, save_path, func, overwrite, skip_if_error, varargin{:});
        WaitMessage.Send;
    end
end
WaitMessage.Destroy;
end

function main_loop(basepath, save_path, func, overwrite, skip_if_error, varargin)
% main_loop: file management & run function
% Inputs:
%    basepath: str path to session
%    save_path: str path to save results
%    func: function to run on each basepath in df (see run)
%    overwrite: bool whether to overwrite existing files in save_path
%    skip_if_error: bool whether to skip if an error occurs
%    varargin: cell array of additional arguments to pass to func

% Get file name from basepath
save_file = encode_file_path(basepath, save_path);

% If file exists and overwrite is False, skip
if exist(save_file, 'file') && ~overwrite
    return;
end

% Calculate some features
if skip_if_error
    try
        results = feval(func, basepath, varargin{:});
    catch
        disp(['Error in ', basepath]);
        return;
    end

else
    results = feval(func, basepath, varargin{:});
end

% Save file
save(save_file, 'results');
end

function save_file = encode_file_path(basepath, save_path)
% encode_file_path: encode file path to be used as a file name
%
% Inputs:
%    basepath: str path to session
%    save_path: str path to save results
%
% Returns:
%    save_file: str encoded file path
%
% example:
%    basepath = 'Z:\Data\AYAold\AB3\AB3_38_41';
%    save_path = 'Z:\home\ryanh\projects\ripple_heterogeneity\replay_02_17_23';
%    batch_analysis.encode_file_path(basepath, save_path)
%    >> 'Z:\home\ryanh\projects\ripple_heterogeneity\replay_02_17_23\Z---___Data___AYAold___AB3___AB3_38_41.mat'

% Normalize paths
basepath = strrep(basepath, filesep, filesep);
save_path = strrep(save_path, filesep, filesep);
% Encode file path with unlikely characters
save_file = fullfile(save_path, [strrep(strrep(basepath, filesep, '___'), ':', '---'), '.mat']);
end

function basepath = decode_file_path(save_file)
% decode_file_path: decode file path to be used as a file name
%
% Inputs:
%    save_file: str encoded file path
%
% Returns:
%    basepath: str path to session
%
% example:
%    save_file = 'Z:\home\ryanh\projects\ripple_heterogeneity\replay_02_17_23\Z---___Data___AYAold___AB3___AB3_38_41.mat'
%    batch_analysis.decode_file_path(save_file)
%    >> 'Z:\Data\AYAold\AB3\AB3_38_41'

% Get basepath from save_file
[~, name, ~] = fileparts(save_file);
basepath = strrep(strrep(name, '___', filesep), '---', ':');
end

function results = load_results(save_path, verbose, add_save_file_name)
% load_results: load results (table) from a .mat file
% Inputs:
%    save_path: str path to save results
%    verbose: bool whether to print progress
%    add_save_file_name: bool whether to add a column with the save file name
%
% Returns:
%    results: table with results

warning('not implemented')
% if ~exist(save_path, 'dir')
%     error(['Folder ', save_path, ' does not exist']);
% end
% 
% sessions = dir(fullfile(save_path, '*.mat'));
% 
% results = table();
% 
% for i = 1:length(sessions)
%     if verbose
%         disp(sessions(i).name);
%     end
%     session_file = fullfile(sessions(i).folder, sessions(i).name);
%     loaded_data = load(session_file);
%     results_ = loaded_data.results;
% 
%     if isempty(results_)
%         continue;
%     end
% 
%     if add_save_file_name
%         results_.save_file_name = repmat({sessions(i).name}, height(results_), 1);
%     end
% 
%     results = [results; results_];
% end
end
