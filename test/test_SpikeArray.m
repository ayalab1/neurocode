function tests = test_SpikeArray
    % Create a test suite for the SpikeArray class 
    tests = functiontests(localfunctions);
end

function testConstructor(testCase)
    % Test SpikeArray constructor
    spikes_cell = {[1, 2, 3], [4, 5, 6], [7, 8, 9]};
    spike_array = SpikeArray(spikes_cell);
    verifyEqual(testCase, spike_array.spikes, [1, 2, 3, 4, 5, 6, 7, 8, 9], 'Error: SpikeArray constructor not working');
    verifyEqual(testCase, spike_array.uid, [1, 1, 1, 2, 2, 2, 3, 3, 3], 'Error: UID not working');
    
    uid = [1, 1, 3, 1, 3];
    spike_array = SpikeArray([1, 2, 3, 4, 5], uid);
    verifyEqual(testCase, spike_array.spikes, [1, 2, 3, 4, 5], 'Error: SpikeArray constructor with UID not working');
    verifyEqual(testCase, spike_array.uid, uid, 'Error: UID not working');
end

function testRestrict(testCase)
    % Test restrict method
    interval_array = IntervalArray([3, 5; 7, 10]);
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    restricted_spike_array = spike_array(interval_array);
    verifyEqual(testCase, restricted_spike_array.spikes, [3, 4, 5, 7, 8, 9], 'Error: restrict method not working');
    verifyEqual(testCase, restricted_spike_array.uid, [1, 2, 2, 3, 3, 3], 'Error: restrict method UID not working');
end

function testNCells(testCase)
    % Test n_cells method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.n_cells(), 3, 'Error: n_cells not working');
end

function testIds(testCase)
    % Test ids method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.ids(), [1, 2, 3]', 'Error: ids not working');
end

function testNSpikes(testCase)
    % Test n_spikes method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.n_spikes(), [3, 3, 3]', 'Error: n_spikes not working');
end

function testFirstEvent(testCase)
    % Test first_event method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.first_event(), 1, 'Error: first_event not working');
end

function testLastEvent(testCase)
    % Test last_event method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.last_event(), 9, 'Error: last_event not working');
end

function testIsSorted(testCase)
    % Test issorted method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.issorted(), true, 'Error: issorted not working');
end

function testIsEmpty(testCase)
    % Test isempty method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    verifyEqual(testCase, spike_array.isempty(), false, 'Error: isempty not working');
end

function testCellSelection(testCase)
    % Test cell selection method
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    first_cell_array = spike_array(1);
    verifyEqual(testCase, first_cell_array.n_cells(), 3, 'Error: cell selection n_cells not working');
    verifyEqual(testCase, first_cell_array.n_active_cells(), 1, 'Error: cell selection n_active_cells not working');
    verifyEqual(testCase, first_cell_array.ids(), [1, 2, 3]', 'Error: cell selection ids not working');
    verifyEqual(testCase, first_cell_array.n_spikes(), [3, 0, 0]', 'Error: cell selection n_spikes not working');
    verifyEqual(testCase, first_cell_array.first_event(), 1, 'Error: cell selection first_event not working');
    verifyEqual(testCase, first_cell_array.last_event(), 3, 'Error: cell selection last_event not working');
    verifyEqual(testCase, first_cell_array.issorted(), true, 'Error: cell selection issorted not working');
    verifyEqual(testCase, first_cell_array.isempty(), false, 'Error: cell selection isempty not working');
end

function testBin(testCase)
    % Test bin
    spike_array = SpikeArray({[1, 2, 3], [4, 5, 6], [7, 8, 9]});
    bst = spike_array.bin('ds', 1);
    X = [1, 0, 0;...
         1, 0, 0;...
         1, 0, 0;...
         0, 1, 0;...
         0, 1, 0;...
         0, 1, 0;...
         0, 0, 1;...
         0, 0, 1;...
         0, 0, 1];
    verifyEqual(testCase, bst.data, X, 'Error: bin method not working');
end

function testToCellArray(testCase)
    % Test to_cell_array
    spikes_cell = {[1, 2, 3], [4, 5, 6], [7, 8, 9]};
    spike_array = SpikeArray(spikes_cell);
    spikes_cell_2 = spike_array.to_cell_array();
    verifyEqual(testCase, spikes_cell, spikes_cell_2, 'Error: to_cell_array not working');
end
