classdef IntervalArray < handle
    % IntervalArray is a class for representing an array of intervals
    %
    % Properties:
    %   intervals - array of intervals represented as [start, stop]
    %
    % Methods:
    %   IntervalArray - constructor, creates a new IntervalArray object
    %   validate_intervals - method that checks if intervals are valid
    %   sort - sorts intervals by start time
    %   starts - returns the start time of each interval
    %   stops - retruns the stop time of each interval
    %   min - minimum bound of all intervals
    %   max - maximum bound of all intervals
    %   is_finite - if min and max are finite
    %   centers - mean (center) val of each interval
    %   n_intervals - returns the number of intervals in the object
    %   expand - expands or shrinks the intervals by a certain amount
    %   isempty - checks if the intervals are empty
    %   lengths - returns the duration of each interval
    %   duration - returns the total duration of all the intervals
    %   intersect - returns the intersection of the current IntervalArray
    %               object and another IntervalArray object
    %   union - returns the union of the current IntervalArray object and
    %           another IntervalArray object
    %   setdiff - returns the set difference of the current IntervalArray
    %             object and another IntervalArray object
    %   complement - returns the complement of the current IntervalArray
    %                object
    %   merge - returns the merger of the current IntervalArray object
    %   plus - returns the addition of the current IntervalArray object
    %          and another IntervalArray object
    %   remove_empty - remove empty intervals
    %   eq - check if intervals are equal
    %   in - check if point is within one of the intervals
    %   plot - plot intervals using PlotIntervals
    %
    % Examples:
    %   myIntervalArray = IntervalArray([0,5;10,15])
    %   myIntervalArray.validate_intervals()
    %   myIntervalArray.sort()
    %   myIntervalArray(1)
    %   myIntervalArray.starts()
    %   myIntervalArray.stops()
    %   myIntervalArray.n_intervals()
    %   myIntervalArray.expand(1)
    %   myIntervalArray.isempty()
    %   myIntervalArray.lengths()
    %   myIntervalArray.duration()
    %   myIntervalArray.intersect(otherIntervalArray)
    %   myIntervalArray.union(otherIntervalArray)
    %   myIntervalArray.setdiff(otherIntervalArray)
    %   myIntervalArray.complement()
    %   myIntervalArray.merger()
    %   myIntervalArray.plus(otherIntervalArray)
    %   myIntervalArray.remove_empty()
    %   myIntervalArray.eq(otherIntervalArray)
    %   myIntervalArray.in(point)
    %   myIntervalArray + otherIntervalArray
    %   myIntervalArray - otherIntervalArray
    %   myIntervalArray & otherIntervalArray
    %   myIntervalArray | otherIntervalArray
    %   ~myIntervalArray

    % Ryan H 2023

    properties
        intervals
    end

    methods
        function obj = IntervalArray(intervals_in)
            if ~exist('intervals_in', 'var')
                intervals_in = [-inf, inf];
            end
            obj.intervals = intervals_in;
            obj.validate_intervals();
            obj.sort();
        end

        function obj = validate_intervals(obj)

            % remove nan
            nan_idx = any(isnan(obj.intervals), 2);
            if any(nan_idx)
                warning('removing intervals with nans')
                obj.intervals(nan_idx, :) = [];
            end

            % check if intervals are valid
            if any(obj.intervals(:, 1) > obj.intervals(:, 2))
                error('Invalid intervals: start time must be less than end time')
            end
        end

        function interval = subsref(obj, S)
            if isequal(S.type, '()')
                if S.subs{1} > obj.n_intervals() || S.subs{1} < 0
                    error('Index out of bounds')
                end
                interval = IntervalArray(obj.intervals(S.subs{1}, :));
            else
                interval = builtin('subsref', obj, S);
            end
        end

        function intervals = data(obj)
            % can access intervals via .data
            intervals = obj.intervals;
        end

        function disp(obj)

            % pull out intervals and remove inf
            intervals_ = obj.intervals;
            intervals_(any(isinf(intervals_), 2), :) = [];

            % calc total duration
            obj_duration = seconds(sum(intervals_(:, 2)-intervals_(:, 1)));

            if obj_duration < seconds(1)
                duration_str = datestr(obj_duration, 'FFF');
                units = 'ms';
            elseif obj_duration < seconds(60)
                duration_str = datestr(obj_duration, 'SS.FFF');
                units = 'seconds';
            elseif obj_duration < seconds(3600)
                duration_str = datestr(obj_duration, 'MM:SS.FFF');
                units = 'minutes';
            elseif obj_duration < seconds(86400)
                duration_str = datestr(obj_duration, 'HH:MM:SS.FFF');
                units = 'hours';
            else
                duration_str = datestr(obj_duration, 'DD:HH:MM:SS.FFF');
                units = 'days';
            end
            fprintf('<%s %d epochs> of length %s %s \n', ...
                "IntervalArray:", ...
                obj.n_intervals, ...
                duration_str, ...
                units);
        end

        function max_ = max(obj)
            % maximum bound of all intervals in IntervalArray
            max_ = max(obj.intervals(:, 2));
        end

        function min_ = min(obj)
            % minimum bound of all intervals in IntervalArray
            min_ = min(obj.intervals(:, 1));
        end

        function is_finite_ = is_finite(obj)
            % Is the interval [start, stop) finite.
            is_finite_ = ~(isinf(obj.min) | isinf(obj.max));
        end

        function centers_ = centers(obj)
            centers_ = mean(obj.intervals, 2);
        end

        function obj = sort(obj)
            % sort intervals by start time
            obj.intervals = sortrows(obj.intervals, 1);
        end

        function new = remove_empty(obj)
            new = IntervalArray();
            new.intervals = obj.intervals;
            % remove empty intervals
            new.intervals(obj.intervals(:, 1) == obj.intervals(:, 2), :) = [];
        end

        function equal = eq(obj, other)
            % check if intervals are equal
            equal = isequal(obj.intervals, other.intervals);
        end

        function contains = in(obj, point)
            % check if point is within one of the intervals
            contains = any(point >= obj.intervals(:, 1) & point <= obj.intervals(:, 2));
        end

        function starts_ = starts(obj)
            % check if intervals are not empty
            if ~obj.isempty()
                starts_ = obj.intervals(:, 1);
            else
                starts_ = [];
            end
        end

        function stops_ = stops(obj)
            % check if intervals are not empty
            if ~obj.isempty()
                stops_ = obj.intervals(:, 2);
            else
                stops_ = [];
            end
        end

        function n_intervals_ = n_intervals(obj)
            % check if intervals are not empty
            if ~obj.isempty()
                n_intervals_ = size(obj.intervals, 1);
            else
                n_intervals_ = 0;
            end
        end

        function new = expand(obj, amount, direction)
            if ~exist('direction', 'var')
                direction = 'both';
            end
            if contains(direction, 'both')
                resize_starts = obj.intervals(:, 1) - amount;
                resize_stops = obj.intervals(:, 2) + amount;
            elseif contains(direction == 'start')
                resize_starts = obj.intervals(:, 1) - amount;
                resize_stops = obj.intervals(:, 2);
            elseif contains(direction == 'stop')
                resize_starts = obj.intervals(:, 1);
                resize_stops = obj.intervals(:, 2) + amount;
            else
                error("direction must be 'both', 'start', or 'stop'")
            end
            new = IntervalArray();
            new.intervals = [resize_starts, resize_stops];
        end

        function isempty_ = isempty(obj)
            % check if intervals are not empty
            isempty_ = isempty(obj.intervals);
        end

        function lengths_ = lengths(obj)
            intervals_ = obj.intervals;
            % check if intervals are not empty
            if ~isempty(intervals_)
                lengths_ = intervals_(:, 2) - intervals_(:, 1);
            else
                lengths_ = [];
            end
        end

        function duration_ = duration(obj)
            % check if intervals are not empty
            if ~isempty(obj.intervals)
                duration_ = sum(obj.lengths);
            else
                duration_ = [];
            end
        end

        function new = and(obj, other)
            % intersection using interval_1 & interval_2
            % https://www.mathworks.com/help/matlab/matlab_oop/implementing-operators-for-your-class.html
            new = intersect(obj, other);
        end

        function new = or(obj, other)
            % union using interval_1 | interval_2
            new = union(obj, other);
        end

        function new = minus(obj, other)
            % union using interval_1 - interval_2
            new = setdiff(obj, other);
        end

        function new = not(obj)
            % complement using ~interval_1
            new = complement(obj);
        end

        function new = intersect(obj, other)
            new = IntervalArray();
            if isa(other, 'IntervalArray')
                % Get intervals
                intervals_ = obj.intervals;
                other_intervals = other.intervals;
                % Initialize variables to store the intersection intervals
                intersection = [];
                i = 1;
                j = 1;
                % Iterate through intervals and other_intervals
                while (i <= size(intervals_, 1) && j <= size(other_intervals, 1))
                    % Check for intersection
                    if (intervals_(i, 2) >= other_intervals(j, 1) && intervals_(i, 1) <= other_intervals(j, 2))
                        start = max(intervals_(i, 1), other_intervals(j, 1));
                        stop = min(intervals_(i, 2), other_intervals(j, 2));
                        intersection = [intersection; start, stop];
                    end
                    % Move to the next interval or other_interval
                    if (intervals_(i, 2) < other_intervals(j, 2))
                        i = i + 1;
                    else
                        j = j + 1;
                    end
                end
                new.intervals = intersection;
            else
                error("unsupported operand type(s) for intersection: %s and %s", ...
                    class(obj), class(other));
            end
        end


        function new = union(obj, other)
            if isa(other, 'IntervalArray')
                intervals_ = obj.intervals;
                other_intervals = other.intervals;
                % check if intervals and other_intervals are not empty
                if ~isempty(intervals_) && ~isempty(other_intervals)
                    i = 1;
                    j = 1;
                    while (i <= size(intervals_, 1) && j <= size(other_intervals, 1))
                        if (intervals_(i, 1) < other_intervals(j, 1))
                            if (intervals_(i, 2) < other_intervals(j, 1))
                                i = i + 1;
                            else
                                if (intervals_(i, 2) < other_intervals(j, 2))
                                    intervals_(i, 2) = other_intervals(j, 2);
                                    j = j + 1;
                                else
                                    j = j + 1;
                                end
                            end
                        else
                            if (other_intervals(j, 2) < intervals_(i, 1))
                                j = j + 1;
                            else
                                if (other_intervals(j, 2) < intervals_(i, 2))
                                    intervals_ = [intervals_(1:i-1, :); ...
                                        other_intervals(j, :); intervals_(i:end, :)];
                                    i = i + 1;
                                    j = j + 1;
                                else
                                    intervals_(i, 1) = other_intervals(j, 1);
                                    j = j + 1;
                                end
                            end
                        end
                    end
                    if (j <= size(other_intervals, 1))
                        intervals_ = [intervals_; other_intervals(j:end, :)];
                    end
                elseif isempty(intervals_)
                    intervals_ = other_intervals;
                end
                new = IntervalArray(intervals_);
            else
                error("unsupported operand type(s) for union: %s and %s", ...
                    class(obj), class(other));
            end
        end

        function new = complement(obj)
            new = IntervalArray();
            % Get intervals
            intervals_ = obj.intervals;
            % Check if intervals are not empty
            if ~isempty(intervals_)
                % Initialize variables to store the complement intervals
                complement = [-inf, intervals_(1, 1)];
                for i = 1:size(intervals_, 1) - 1
                    complement = [complement; intervals_(i, 2), intervals_(i+1, 1)];
                end
                complement = [complement; intervals_(end, 2), inf];
                new.intervals = complement;
            else
                new.intervals = [-inf, inf];
            end
        end

        function new = merge(obj, varargin)
            new = IntervalArray();
            if isempty(varargin)
                [new.intervals, ~] = ConsolidateIntervals(obj.intervals);
            else
                [new.intervals, ~] = ConsolidateIntervals(obj.intervals, varargin);
            end
        end

        function new = plus(obj, other)
            if isa(other, 'IntervalArray')
                new = IntervalArray();
                intervals_ = obj.intervals;
                other_intervals = other.intervals;
                % check if intervals and other_intervals are not empty
                if ~isempty(intervals_) && ~isempty(other_intervals)
                    new_intervals = [intervals_; other_intervals];
                    new_intervals = sortrows(new_intervals);
                    new.intervals = new_intervals;
                elseif isempty(intervals_)
                    new = other;
                end
            else
                error("unsupported operand type(s) for +: %s and %s", ...
                    class(obj), class(other));
            end
        end

        function new = setdiff(obj, other)
            if isa(other, 'IntervalArray')
                new = IntervalArray();
                other_intervals = other.intervals;
                intervals_ = obj.intervals;
                % loop through the intervals of other
                for i = 1:size(other_intervals, 1)
                    % loop through the intervals of the object
                    j = 1;
                    for jj = 1:size(intervals_, 1)
                        if (intervals_(j, 1) >= other_intervals(i, 1) && ...
                                intervals_(j, 2) <= other_intervals(i, 2))
                            % interval is completely inside the other interval
                            intervals_(j, :) = [];
                            continue
                        elseif (intervals_(j, 1) < other_intervals(i, 1) && ...
                                intervals_(j, 2) > other_intervals(i, 1) && ...
                                intervals_(j, 2) <= other_intervals(i, 2))
                            % interval starts before and ends inside the other interval
                            intervals_(j, 2) = other_intervals(i, 1);
                        elseif (intervals_(j, 1) >= other_intervals(i, 1) && ...
                                intervals_(j, 1) < other_intervals(i, 2) && ...
                                intervals_(j, 2) > other_intervals(i, 2))
                            % interval starts inside and ends after the other interval
                            intervals_(j, 1) = other_intervals(i, 2);
                        elseif (intervals_(j, 1) < other_intervals(i, 1) && ...
                                intervals_(j, 2) > other_intervals(i, 2))
                            % interval starts before and ends after the other interval
                            new_interval = [other_intervals(i, 2), intervals_(j, 2)];
                            intervals_(j, 2) = other_intervals(i, 1);
                            intervals_ = [intervals_; new_interval];
                        end
                        j = j + 1;
                    end
                end
                new.intervals = intervals_;
            else
                error("unsupported operand type(s) for setdiff: %s and %s", ...
                    class(obj), class(other));
            end
        end

        function out = plot(obj, varargin)
            if isempty(varargin)
                yLim = ylim;
                alphaValue = 0.5;
                dy = yLim(2) - yLim(1);
                colors = rand(size(obj.intervals, 1), 3);
                for i = 1:obj.n_intervals
                    % Better off implementing "PlotIntervals" directly here to avoid calling "uistack" multiple times
                    dx = diff(obj.intervals(i, :));
                    out(i) = patch(obj.intervals(i, 1)+[0, 0, dx, dx], yLim(1)+[0, dy, dy, 0], colors(i, :), 'LineStyle', 'none');
                    alpha(out(i), alphaValue);
                end
                uistack(out, 'bottom');
            else
                out = PlotIntervals(obj.intervals, varargin);
            end
        end
    end
end