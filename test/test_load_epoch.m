function tests = test_load_epoch
tests = functiontests(localfunctions);
end

function test_basic_epoch_loading(testCase)

% Create temporary directory
basepath = tempname;
mkdir(basepath);

% Use real basename convention
[~, basename] = fileparts(basepath);

% Create fake session struct
session.epochs = {; ...
    struct( ...
    'name', 'sleep', ...
    'startTime', int32(0), ...
    'stopTime', int32(100), ...
    'environment', 'home', ...
    'behavioralParadigm', 'none', ...
    'manipulation', 'none', ...
    'stimuli', 'none', ...
    'notes', 'baseline' ...
    ), ...
    struct( ...
    'name', 'run', ...
    'startTime', single(100), ...
    'stopTime', single(200) ...
    ); ...
    };

% Save session file with correct basename
save(fullfile(basepath, [basename, '.session.mat']), 'session');

% Run function
epochs = load_epoch('basepath', basepath);

% ---- Assertions ----

verifyEqual(testCase, height(epochs), 2);

verifyEqual(testCase, epochs.name{1}, 'sleep');
verifyEqual(testCase, epochs.name{2}, 'run');

verifyClass(testCase, epochs.startTime, 'double');
verifyClass(testCase, epochs.stopTime, 'double');

verifyEqual(testCase, epochs.startTime(1), 0);
verifyEqual(testCase, epochs.stopTime(2), 200);

verifyEqual(testCase, epochs.environment{1}, 'home');
verifyEqual(testCase, epochs.environment{2}, '');

verifyEqual(testCase, epochs.behavioralParadigm{2}, '');
verifyEqual(testCase, epochs.manipulation{2}, '');
verifyEqual(testCase, epochs.stimuli{2}, '');
verifyEqual(testCase, epochs.notes{2}, '');

verifyEqual(testCase, epochs.basepath{1}, basepath);
verifyEqual(testCase, epochs.basepath{2}, basepath);

end
