function test_regress_out()
    % Test case 1
    a = [1 2 3 4 5];
    b = [2 4 6 8 10];
    expected_result = double([3.0 3.0 3.0 3.0 3.0]);
    result = regress_out(a, b);
    assert(all(abs(result-expected_result) < 1e4*eps(min(abs(result),abs(expected_result)))), 'Test case 1 failed.');

    % Test case 2
    a = [1 2 3 4 5];
    b = [0 0 0 0 0];
    expected_result = [1 2 3 4 5];
    result = regress_out(a, b);
    assert(all(abs(result-expected_result) < 1e4*eps(min(abs(result),abs(expected_result)))), 'Test case 2 failed.');

    % Test case 3
    a = [1 2 3 4 5];
    b = [1 1 1 1 1];
    expected_result = [1 2 3 4 5];
    result = regress_out(a, b);
    assert(all(abs(result-expected_result) < 1e4*eps(min(abs(result),abs(expected_result)))), 'Test case 3 failed.');

    % Add more test cases as needed

    disp('All test cases passed.');
end