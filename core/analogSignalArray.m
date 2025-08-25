classdef analogSignalArray < handle
    % analogSignalArray is a class for representing continous single or
    %   multi-channel signals
    %
    % analogSignalArray can be position tracking, lfp, or any continous
    % signal
    %
    % Properties:
    %   data - array of data represented as [n samples x n signals]
    %   timestamps - array of timestamps represented as [n samples]
    %   sampling_rate - sample rate of data (will estimate if not provided)
    %
    % Methods:
    %   analogSignalArray - constructor, creates a new analogSignalArray object
    %   validate_signals - method that checks if data is valid
    %   restrict - restrict data to the times within an IntervalArray object
    %   duration - returns duration of signal
    %   n_signals - returns number of signals
    %   n_samples - returns number of samples
    %   isempty - checks if the analogSignalArray is empty
    %   issorted - checks if the analogSignalArray is sorted (time)
    %   smooth - smooth signal with gaussian kernel
    %
    % Example:
    %
    % First, lets put session epochs into IntervalArray
    %
    % start = [];
    % stop = [];
    % for ep = 1:length(session.epochs)
    %     start = [start,session.epochs{ep}.startTime];
    %     stop = [stop,session.epochs{ep}.stopTime];
    % end
    %
    % Now, put tracking data into analogSignalArray
    %
    % positions = analogSignalArray(...
    %     'data',[behavior.position.x;behavior.position.y],...
    %     'timestamps',behavior.timestamps)
    %
    % Lets restrict to epoch 3
    %
    % positions(epochs(3))
    % ans =
    % <analogSignalArray: 2 signals> of length 54:06.362 minutes
    %
    % check if same amount of time in epoch 3
    %
    %  epochs(3)
    % <IntervalArray: 1 epochs> of length 54:06.387 minutes
    %
    % Ryan Harvey 2023

    properties
        data
        timestamps
        sampling_rate
    end

    methods
        function self = analogSignalArray(varargin)
            self.data = [];
            self.timestamps = [];

            if nargin == 0
                return;
            end

            p = inputParser;
            addParameter(p, 'data', []);
            addParameter(p, 'timestamps', []);
            addParameter(p, 'sampling_rate', []);

            parse(p, varargin{:});
            self.data = p.Results.data;
            self.timestamps = p.Results.timestamps;
            self.sampling_rate = p.Results.sampling_rate;

            self = validate_signals(self);
        end

        function self = validate_signals(self)

            % make into (sample x signal) format
            n_rows_data = size(self.data, 1);
            n_cols_data = size(self.data, 2);

            n_rows_timestamps = size(self.timestamps, 1);
            n_cols_timestamps = size(self.timestamps, 2);

            n_timestamps = max(n_rows_timestamps, n_cols_timestamps);

            % [n samples x n signals]
            if n_cols_data == n_timestamps
                self.data = self.data';
            elseif n_rows_data == n_timestamps
                % great, correct shape
            else
                error('timestamps and data have different n samples')
            end

            if n_rows_timestamps < n_cols_timestamps
                self.timestamps = self.timestamps';
            end

            % check if data and timestamps are the same length
            n_rows_data = size(self.data, 1);
            n_rows_timestamps = size(self.timestamps, 1);

            if n_rows_data ~= n_rows_timestamps
                error('timestamps and data have different n samples')
            end

            % make sure is sorted
            if ~self.issorted()
                [self.timestamps, sort_idx] = sort(self.timestamps);
                self.data = self.data(sort_idx, :);
            end

            % if no sample rate is provided, estimate from timestamps
            if isempty(self.sampling_rate)
                self.sampling_rate = 1 / mode(diff(self.timestamps));
            end
        end

        function asa = subsref(self, S)
            if isequal(S.type, '()')
                if isa(S.subs{1}, 'IntervalArray')
                    asa = restrict(self, S.subs{1});
                end
            else
                asa = builtin('subsref', self, S);
            end
        end

        function out = fs(self)
            out = self.sampling_rate;
        end

        function asa = restrict(self, intervals, varargin)
            asa = analogSignalArray();
            if isempty(varargin)
                [asa.timestamps, idx] = Restrict(self.timestamps, intervals.intervals);
            else
                [asa.timestamps, idx] = Restrict(self.timestamps, intervals.intervals, varargin{:});
            end
            asa.data = self.data(idx, :);
            asa.sampling_rate = self.sampling_rate;
        end

        function duration_ = duration(self)
            % calculate duration of contiguous timestamps
            ts_diff = diff(self.timestamps);
            ts_diff(ts_diff > (1 / self.sampling_rate)*2) = [];
            duration_ = sum(ts_diff);
        end

        function disp(self)
            obj_duration = seconds(self.duration);
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

            fprintf('<%s %d signals %.3f Hz> of length %s %s \n', ...
                "analogSignalArray:", ...
                self.n_signals, ...
                self.sampling_rate, ...
                duration_str, ...
                units);
        end

        function n_signals_ = n_signals(self)
            n_signals_ = size(self.data, 2);
        end

        function n_samples_ = n_samples(self)
            n_samples_ = size(self.data, 1);
        end

        function isempty_ = isempty(self)
            isempty_ = isempty(self.data);
        end

        function sorted = issorted(self)
            sorted = all(sort(self.timestamps) == sort(self.timestamps));
        end

        function asa = smooth(self, varargin)
            % gaussian smooth using Smooth.m
            p = inputParser;
            addParameter(p, 'window', 0.05); % standard deviations (50ms default)
            parse(p, varargin{:});
            window = p.Results.window;

            % convert window into samples
            window = window * self.sampling_rate;

            % remove window from varargin to avoid error with Smooth
            idx = contains(varargin{:, 1}, 'window');
            varargin(idx, :) = [];

            % call Smooth and pass args
            asa.data = Smooth(self.data, ...
                [window, 0], ...
                varargin{:});
        end

    end
end
