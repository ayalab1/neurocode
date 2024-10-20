function tests = test_regress_out
% Create a test suite for the regress_out function
tests = functiontests(localfunctions);
end

function testCase1(testCase)
a = [1, 2, 3, 4, 5];
b = [2, 4, 6, 8, 10];
expected_result = double([3.0, 3.0, 3.0, 3.0, 3.0]);
result = regress_out(a, b);
verifyEqual(testCase, result, expected_result, 'AbsTol', 1e-4, 'Test case 1 failed.');
end

function testCase2(testCase)
a = [1, 2, 3, 4, 5];
b = [0, 0, 0, 0, 0];
expected_result = [1, 2, 3, 4, 5];
result = regress_out(a, b);
verifyEqual(testCase, result, expected_result, 'AbsTol', 1e-4, 'Test case 2 failed.');
end

function testCase3(testCase)
a = [1, 2, 3, 4, 5];
b = [1, 1, 1, 1, 1];
expected_result = [1, 2, 3, 4, 5];
result = regress_out(a, b);
verifyEqual(testCase, result, expected_result, 'AbsTol', 1e-4, 'Test case 3 failed.');
end

% Add more test functions as needed
