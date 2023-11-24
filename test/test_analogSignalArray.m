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
