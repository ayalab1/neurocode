function tests = test_analogSignalArray
    tests = functiontests(localfunctions);
end

function testConstructor(testCase)
    data = [1 2 3; 4 5 6];
    timestamps = [0 1 2];
    sampling_rate = 1;

    asa = analogSignalArray('data', data, 'timestamps', timestamps,...
        'sampling_rate', sampling_rate);

    verifyEqual(testCase, asa.data, data');
    verifyEqual(testCase, asa.timestamps, timestamps');
    verifyEqual(testCase, asa.sampling_rate, sampling_rate);
end

function testRestrict(testCase)
    data = [1 2 3 4 5 6 7 8 9 10 11; 1 2 3 4 5 6 7 8 9 10 11];
    timestamps = [0 1 2 3 4 5 6 7 8 9 10];
    intervalArray = IntervalArray([0,2 ; 4,6]);
    asa = analogSignalArray('data', data, 'timestamps', timestamps);

    restricted_asa = asa.restrict(intervalArray);
    
    timestamps = [0,1,2,4,5,6]';
    data = [1,2,3,5,6,7;1,2,3,5,6,7]';
    verifyEqual(testCase, restricted_asa.data, data);
    verifyEqual(testCase, restricted_asa.timestamps, timestamps);
end

function testDuration(testCase)
    data = [1 2 3; 4 5 6];
    timestamps = [0 1 2];
    asa = analogSignalArray('data', data, 'timestamps', timestamps);

    verifyEqual(testCase, asa.duration(), 2);
end

function testN_Signals(testCase)
    data = [1 2 3; 4 5 6];
    timestamps = [0 1 2];
    asa = analogSignalArray('data', data, 'timestamps', timestamps);

    verifyEqual(testCase, asa.n_signals(), 2);
end

function testN_Samples(testCase)
    data = [1 2 3; 4 5 6];
    timestamps = [0 1 2];
    asa = analogSignalArray('data', data, 'timestamps', timestamps);

    verifyEqual(testCase, asa.n_samples(), 3);
end

function testIsempty(testCase)
    data = [1 2 3; 4 5 6];
    timestamps = [0 1 2];
    asa = analogSignalArray('data', data, 'timestamps', timestamps);

    verifyEqual(testCase, asa.isempty(), false);
end

function testIssorted(testCase)
    data = [1 2 3; 4 5 6];
    timestamps = [0 1 2];
    asa = analogSignalArray('data', data, 'timestamps', timestamps);

    verifyEqual(testCase, asa.issorted(), true);
end

function testDownsample(testCase)
    % Create synthetic signal
    fs = 100; % original sampling rate
    t = (0:1/fs:5)'; % 5 seconds
    data = [sin(2*pi*5*t), cos(2*pi*2*t)]; % 2-channel signal

    asa = analogSignalArray( ...
        'data', data, ...
        'timestamps', t, ...
        'sampling_rate', fs);

    % Downsample
    new_fs = 10;
    asa.downsample(new_fs);

    % Check sampling rate updated
    verifyEqual(testCase, asa.sampling_rate, new_fs);

    % Check timestamps spacing
    dt = diff(asa.timestamps);
    verifyEqual(testCase, mean(dt), 1/new_fs, 'AbsTol', 1e-10);

    % Check number of channels preserved
    verifyEqual(testCase, size(asa.data, 2), 2);

    % Check number of samples roughly correct
    expected_n = floor(length(t) * (new_fs / fs));
    verifyLessThanOrEqual(testCase, abs(size(asa.data,1) - expected_n), 1);

    % Check timestamps start correctly
    verifyEqual(testCase, asa.timestamps(1), t(1));
end

