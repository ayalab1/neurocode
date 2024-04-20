function test_IntervalArray()
    % Test constructor
    myIntervalArray = IntervalArray([0,5;10,15]);
    assert(isequal(myIntervalArray.intervals,[0,5;10,15]),'Error: constructor not working');

    % Test validate_intervals
    myIntervalArray.intervals = [10,5;10,15];
    try
        myIntervalArray.validate_intervals();
        error('Error: validate_intervals not working')
    catch
    end

    % Test sort
    myIntervalArray.intervals = [10,15;0,5];
    myIntervalArray.sort();
    assert(isequal(myIntervalArray.intervals,[0,5;10,15]),'Error: sort not working');

    % Test starts
    assert(isequal(myIntervalArray.starts(),[0;10]),'Error: starts not working');

    % Test stops
    assert(isequal(myIntervalArray.stops(),[5;15]),'Error: stops not working');

    % Test n_intervals
    assert(isequal(myIntervalArray.n_intervals(),2),'Error: n_intervals not working');

    % Test expand
    myIntervalArray = myIntervalArray.expand(1);
    assert(isequal(myIntervalArray.intervals,[-1,6;9,16]),'Error: expand not working');

    % Test isempty
    assert(isequal(myIntervalArray.isempty(),false),'Error: isempty not working');

    % Test lengths
    assert(isequal(myIntervalArray.lengths(),[7;7]),'Error: lengths not working');

    % Test duration
    assert(isequal(myIntervalArray.duration(),14),'Error: duration not working');

    % Test intersect
    otherIntervalArray = IntervalArray([-2,2;8,12]);
    myIntervalArray = myIntervalArray.intersect(otherIntervalArray);
    assert(isequal(myIntervalArray.intervals,[-1,2;9,12]),'Error: intersect not working');

    % Test union
    otherIntervalArray = IntervalArray([-2,2;8,12]);
    myIntervalArray = myIntervalArray.union(otherIntervalArray);
    assert(isequal(myIntervalArray.intervals,[-2,2;8,12]),'Error: union not working');

    % Test setdiff
    otherIntervalArray = IntervalArray([-2,2;18,20]);
    myIntervalArray = myIntervalArray.setdiff(otherIntervalArray);
    assert(isequal(myIntervalArray.intervals,[8,12]),'Error: setdiff not working');

    % Test complement
    myIntervalArray = myIntervalArray.complement();
    assert(isequal(myIntervalArray.intervals,[-inf,8;12,inf]),'Error: complement not working');

    % Test merge
    myIntervalArray = IntervalArray([10,20;15,25]);
    myIntervalArray = myIntervalArray.merge();
    assert(isequal(myIntervalArray.intervals,[10,25]),'Error: merge not working');

    % Test plus
    otherIntervalArray = IntervalArray([-2,2;8,12]);
    myIntervalArray = myIntervalArray.plus(otherIntervalArray);
    assert(isequal(myIntervalArray.intervals,[-2,2;8,12;10,25]),'Error: plus not working');
    
    % Test remove_empty
    myIntervalArray = IntervalArray([0,2;5,5;8,12]);
    myIntervalArray = myIntervalArray.remove_empty();
    assert(isequal(myIntervalArray.intervals,[0,2;8,12]),'Error: remove_empty not working');
    
    % Test eq
    myIntervalArray = IntervalArray([0,2;8,12]);
    otherIntervalArray = IntervalArray([0,2;8,12]);
    is_equal = myIntervalArray.eq(otherIntervalArray);
    assert(isequal(is_equal,true),'Error: eq not working');
    
    % Test in
    myIntervalArray = IntervalArray([0,2;8,12]);
    point = 9;
    contain = myIntervalArray.in(point);
    assert(isequal(contain,true),'Error: in not working');

    % Test empty interval
    myIntervalArray = IntervalArray([]);
    assert(isequal(myIntervalArray.isempty,true),'Error: empty intervals not working');
end