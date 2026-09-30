function tests = test_ConsolidateIntervals
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.originalPath = path;
addpath(fullfile(root, 'utilities', 'intervalsC+'), ...
    fullfile(root, 'utilities', 'Helpers'), fullfile(root, 'lfp'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testSingleInterval(testCase)
interval = [7426.8, 7433.8];
for strict = {'off', 'on'}
    [consolidated, target] = ConsolidateIntervals(interval, 'strict', strict{1});
    verifyEqual(testCase, consolidated, interval);
    verifyEqual(testCase, target, 1);
end
end

function testSingleColumnInterval(testCase)
interval = [7426.8; 7433.8];
[consolidated, target] = ConsolidateIntervals(interval);
verifyEqual(testCase, consolidated, interval');
verifyEqual(testCase, target, 1);
end

function testEmptyIntervals(testCase)
[consolidated, target] = ConsolidateIntervals(zeros(0, 2));
verifyEqual(testCase, consolidated, zeros(0, 2));
verifyEmpty(testCase, target);
end

function testOverlappingIntervals(testCase)
intervals = [10, 20; 1, 3; 2, 5; 15, 25];
for strict = {'off', 'on'}
    [consolidated, target] = ConsolidateIntervals(intervals, 'strict', strict{1});
    verifyEqual(testCase, consolidated, [1, 5; 10, 25]);
    verifyEqual(testCase, target, [2; 1; 1; 2]);
end
end

function testTouchingIntervals(testCase)
intervals = [1, 2; 2, 3];
[consolidated, target] = ConsolidateIntervals(intervals);
verifyEqual(testCase, consolidated, [1, 3]);
verifyEqual(testCase, target, [1; 1]);
[consolidated, target] = ConsolidateIntervals(intervals, 'strict', 'on');
verifyEqual(testCase, consolidated, intervals);
verifyEqual(testCase, target, [1; 2]);
end

function testCleanLFPSingleArtefact(testCase)
t = (0:0.01:10)';
values = zeros(size(t));
values(t >= 5 & t <= 5.02) = 100;
[clean, bad, badIntervals] = CleanLFP([t, values], 'thresholds', [8, Inf]);
expectedInterval = [4.5, 5.52];
verifyEqual(testCase, badIntervals, expectedInterval, 'AbsTol', 1e-12);
verifyEqual(testCase, bad, t >= expectedInterval(1) & t <= expectedInterval(2));
verifyEqual(testCase, clean, [t, zeros(size(t))]);
end

function testCleanLFPSingleDerivativeArtefact(testCase)
t = (0:0.01:10)';
values = zeros(size(t));
values(t >= 5) = 1;
[clean, bad, badIntervals] = CleanLFP([t, values], 'thresholds', [Inf, 1]);
expectedInterval = [4.89, 5.09];
verifyEqual(testCase, badIntervals, expectedInterval, 'AbsTol', 1e-12);
verifyTrue(testCase, any(bad));
verifyEqual(testCase, clean(~bad, :), [t(~bad), values(~bad)]);
verifyTrue(testCase, all(isfinite(clean(:))));
end

function testCleanLFPWithoutArtefacts(testCase)
t = (0:0.01:10)';
lfp = [t, sin(t)];
[clean, bad, badIntervals] = CleanLFP(lfp, 'thresholds', [Inf, Inf]);
verifyEqual(testCase, clean, lfp);
verifyFalse(testCase, any(bad));
verifyEqual(testCase, badIntervals, zeros(0, 2));
end
