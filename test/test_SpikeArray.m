function test_SpikeArray()
    % Test SpikeArray constructor
    spikes_cell = {[1, 2, 3], [4, 5, 6], [7, 8, 9]};
    spike_array = SpikeArray(spikes_cell);
    assert(isequal(spike_array.spikes, [1, 2, 3, 4, 5, 6, 7, 8, 9]));
    assert(isequal(spike_array.uid, [1, 1, 1, 2, 2, 2, 3, 3, 3]));
    
    uid = [1, 1, 3, 1, 3];
    spike_array = SpikeArray([1, 2, 3, 4, 5], uid);
    assert(isequal(spike_array.spikes, [1, 2, 3, 4, 5]));
    assert(isequal(spike_array.uid, uid));

    % Test restrict method
    spikes_cell = {[1, 2, 3], [4, 5, 6], [7, 8, 9]};
    spike_array = SpikeArray(spikes_cell);
    interval_array = IntervalArray([3, 5; 7, 10]);
    restricted_spike_array = spike_array(interval_array);
    assert(isequal(restricted_spike_array.spikes, [3, 4, 5, 7, 8, 9]));
    assert(isequal(restricted_spike_array.uid, [1, 2, 2, 3, 3, 3]));

    % Test n_cells method
    assert(spike_array.n_cells() == 3);

    % Test ids method
    assert(isequal(spike_array.ids(), [1, 2, 3]'));

    % Test n_spikes method
    assert(isequal(spike_array.n_spikes(), [3, 3, 3]'));

    % Test first_event method
    assert(spike_array.first_event() == 1);

    % Test last_event method
    assert(spike_array.last_event() == 9);

    % Test issorted method
    assert(spike_array.issorted() == 1);

    % Test isempty method
    assert(spike_array.isempty() == 0);
    
    % Test cell selection method
    first_cell_array = spike_array(1);
    assert(first_cell_array.n_cells() == 3);
    assert(first_cell_array.n_active_cells() == 1);
    assert(all(first_cell_array.ids() == [1,2,3]'));
    assert(all(first_cell_array.n_spikes() == [3,0,0]'));
    assert(first_cell_array.first_event() == 1);
    assert(first_cell_array.last_event() == 3);
    assert(first_cell_array.issorted() == 1);
    assert(first_cell_array.isempty() == 0);
    
    % Test bin
    spikes_cell = {[1, 2, 3], [4, 5, 6], [7, 8, 9]};
    spike_array = SpikeArray(spikes_cell);
    bst = spike_array.bin('ds',1);
    X = [1,0,0;...
        1,0,0;...
        1,0,0;...
        0,1,0;...
        0,1,0;...
        0,1,0;...
        0,0,1;...
        0,0,1;...
        0,0,1];
    assert(isequal(bst.data,X))
    
    % Test to_cell_array
    spikes_cell = {[1, 2, 3], [4, 5, 6], [7, 8, 9]};
    spike_array = SpikeArray(spikes_cell);
    spikes_cell_2 = spike_array.to_cell_array();
    assert(isequal(spikes_cell,spikes_cell_2))
end
