function tests = test_IntervalArray
% Create a test suite for the IntervalArray class
tests = functiontests(localfunctions);
end

function testConstructor(testCase)
% Test constructor
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.intervals, [0, 5; 10, 15], 'Error: constructor not working');
end

function testValidateIntervals(testCase)
% Test validate_intervals
myIntervalArray = IntervalArray([0, 5; 10, 15]);
myIntervalArray.intervals = [10, 5; 10, 15];
% Try to validate intervals and catch any error
try
    myIntervalArray.validate_intervals();
    % If no error is thrown, the test should fail
    testCase.verifyFail('Expected an error but did not receive one.');
catch ME
    % Check if the caught error message is the expected one
    testCase.verifyEqual(ME.message, 'Invalid intervals: start time must be less than end time', ...
        'Error: validate_intervals did not throw the expected message.');
end
end

function testSort(testCase)
% Test sort
myIntervalArray = IntervalArray([10, 15; 0, 5]);
myIntervalArray.sort();
verifyEqual(testCase, myIntervalArray.intervals, [0, 5; 10, 15], 'Error: sort not working');
end

function testStarts(testCase)
% Test starts
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.starts(), [0; 10], 'Error: starts not working');
end

function testStops(testCase)
% Test stops
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.stops(), [5; 15], 'Error: stops not working');
end

function testNIntervals(testCase)
% Test n_intervals
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.n_intervals(), 2, 'Error: n_intervals not working');
end

function testExpand(testCase)
% Test expand
myIntervalArray = IntervalArray([0, 5; 10, 15]);
myIntervalArray = myIntervalArray.expand(1);
verifyEqual(testCase, myIntervalArray.intervals, [-1, 6; 9, 16], 'Error: expand not working');
end

function testIsEmpty(testCase)
% Test isempty
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.isempty(), false, 'Error: isempty not working');
end

function testLengths(testCase)
% Test lengths
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.lengths(), [5; 5], 'Error: lengths not working');
end

function testDuration(testCase)
% Test duration
myIntervalArray = IntervalArray([0, 5; 10, 15]);
verifyEqual(testCase, myIntervalArray.duration(), 10, 'Error: duration not working');
end

function testIntersect(testCase)
% Test intersect
myIntervalArray = IntervalArray([0, 5; 10, 15]);
otherIntervalArray = IntervalArray([-2, 2; 8, 12]);
myIntervalArray = myIntervalArray.intersect(otherIntervalArray);
verifyEqual(testCase, myIntervalArray.intervals, [0, 2; 10, 12], 'Error: intersect not working');
end

function testUnion(testCase)
% Test union
myIntervalArray = IntervalArray([0, 5; 10, 15]);
otherIntervalArray = IntervalArray([-2, 2; 8, 12]);
myIntervalArray = myIntervalArray.union(otherIntervalArray);
verifyEqual(testCase, myIntervalArray.intervals, [-2, 5; 8, 15], 'Error: union not working');
end

function testSetDiff(testCase)
% Test setdiff
myIntervalArray = IntervalArray([0, 5; 10, 15]);
otherIntervalArray = IntervalArray([-2, 2; 18, 20]);
myIntervalArray = myIntervalArray.setdiff(otherIntervalArray);
verifyEqual(testCase, myIntervalArray.intervals, [10, 15], 'Error: setdiff not working');
end

function testComplement(testCase)
% Test complement
myIntervalArray = IntervalArray([0, 5; 10, 15]);
myIntervalArray = myIntervalArray.complement();
verifyEqual(testCase, myIntervalArray.intervals, [-inf, 0; 5, 10; 15, inf], 'Error: complement not working');
end

function testMerge(testCase)
% Test merge
myIntervalArray = IntervalArray([10, 20; 15, 25]);
myIntervalArray = myIntervalArray.merge();
verifyEqual(testCase, myIntervalArray.intervals, [10, 25], 'Error: merge not working');
end

function testPlus(testCase)
% Test plus
myIntervalArray = IntervalArray([10, 20; 15, 25]);
otherIntervalArray = IntervalArray([-2, 2; 8, 12]);
myIntervalArray = myIntervalArray.plus(otherIntervalArray);
verifyEqual(testCase, myIntervalArray.intervals, [-2, 2; 8, 25], 'Error: plus not working');
end

function testRemoveEmpty(testCase)
% Test remove_empty
myIntervalArray = IntervalArray([0, 2; 5, 5; 8, 12]);
myIntervalArray = myIntervalArray.remove_empty();
verifyEqual(testCase, myIntervalArray.intervals, [0, 2; 8, 12], 'Error: remove_empty not working');
end

function testEq(testCase)
% Test eq
myIntervalArray = IntervalArray([0, 2; 8, 12]);
otherIntervalArray = IntervalArray([0, 2; 8, 12]);
is_equal = myIntervalArray.eq(otherIntervalArray);
verifyEqual(testCase, is_equal, true, 'Error: eq not working');
end

function testIn(testCase)
% Test in
myIntervalArray = IntervalArray([0, 2; 8, 12]);
point = 9;
contain = myIntervalArray.in(point);
verifyEqual(testCase, contain, true, 'Error: in not working');
end

function testEmptyInterval(testCase)
% Test empty interval
myIntervalArray = IntervalArray([]);
verifyEqual(testCase, myIntervalArray.isempty(), true, 'Error: empty intervals not working');
end
