classdef SpikeArray < handle
    % SpikeArray is a class for representing an array of spike timestamps
    % from multiple neurons. It can be sliced by IntervalArray and can
    % contain timestamps from multiple neurons.
    %
    % Properties:
    %   spikes - array of spike timestamps
    %   uid - unique identifier for each neuron
    %
    % Methods:
    %   SpikeArray - constructor, creates a new SpikeArray object
    %   restrict - restrict spikes to the times within an IntervalArray object
    %   n_cells - returns the number of neurons in the SpikeArray
    %   ids - returns the unique identifiers of each neuron
    %   n_spikes - returns the number of spikes for each neuron
    %   first_event - returns the timestamp of the first spike
    %   last_event - returns the timestamp of the last spike
    %   issorted - checks if the spikes are sorted in chronological order
    %   isempty - checks if the SpikeArray is empty
    %   plot - plots the spikes using RasterPlot
    %
    % Examples:
    %   mySpikeArray = SpikeArray({[1,2,3], [4,5,6], [7,8,9]})
    %   myIntervalArray = IntervalArray([1,5;7,9])
    %   mySpikeArray(myIntervalArray) or mySpikeArray.restrict(myIntervalArray)
    %   mySpikeArray(1)
    %   mySpikeArray.n_cells()
    %   mySpikeArray.ids()
    %   mySpikeArray.n_spikes()
    %   mySpikeArray.first_event()
    %   mySpikeArray.last_event()
    %   mySpikeArray.issorted()
    %   mySpikeArray.isempty()
    %   mySpikeArray.plot()
    %
    % Ryan H 2023
    
    properties
        spikes
        uid
    end
    
    methods
        function obj = SpikeArray(varargin)
            obj.spikes = [];
            obj.uid = [];
            
            if nargin == 0
                return;
            end
            
            if nargin == 1
                spikes_in = varargin{1};
                if isnumeric(spikes_in)
                    obj.spikes = spikes_in;
                    obj.uid = ones(size(spikes_in));
                elseif iscell(spikes_in)
                    numcells = length(spikes_in);
                    for cc = 1:numcells
                        groups{cc} = cc * ones(size(spikes_in{cc}));
                    end
                    if numcells > 0
                        alltimes = cat(1,spikes_in{:});
                        groups_ = cat(1,groups{:});
                        % check if dim of cat was correct and fix if not
                        n_spikes_total = sum(cellfun('length',spikes_in));
                        if n_spikes_total ~= length(alltimes)
                            alltimes = cat(2,spikes_in{:});
                            groups_ = cat(2,groups{:});
                        end
                        [obj.spikes, sortidx] = sort(alltimes);
                        obj.uid = groups_(sortidx);
                    end
                else
                    error('Input must be numeric array or cell array of numeric arrays.');
                end
            elseif nargin == 2
                spikes_in = varargin{1};
                uid_in = varargin{2};
                if isnumeric(spikes_in) && isnumeric(uid_in) && (size(spikes_in, 1) == size(uid_in, 1))
                    [obj.spikes, sortidx] = sort(spikes_in);
                    obj.uid = uid_in(sortidx);
                else
                    error('Input spikes and uid must be numeric arrays with matching number of rows.');
                end
            else
                error('Invalid number of input arguments.');
            end
        end
        
        function st = subsref(obj,S)
            if isequal(S.type,'()')
                if isa(S.subs{1},'IntervalArray')
                    st = restrict(obj, S.subs{1});
                    return
                end
                if ~any(ismember(S.subs{1},obj.ids))
                    error('Index out of bounds')
                end
                st = SpikeArray();
                idx = obj.uid == S.subs{1};
                st.spikes = obj.spikes(idx);
                st.uid = obj.uid(idx);
            else
                st = builtin('subsref',obj,S);
            end
        end

        function st  = restrict(obj, intervals, varargin)
            st  = SpikeArray();
            if isempty(varargin)
                [st.spikes, idx] = Restrict(obj.spikes, intervals.intervals);
            else
                [st.spikes, idx] = Restrict(obj.spikes, intervals.intervals, varargin{:});
            end
            st.uid = obj.uid(idx);
        end
        
        function bst = bin(obj,varargin)
            
            p = inputParser;
            addParameter(p,'ds',0.0625,@isnumeric)
            parse(p,varargin{:})
            ds = p.Results.ds;
            
            % set up bin edge or each cell
            uid_edge = [obj.ids;max(obj.ids)+1] - .5;
            
            [bst,~,~] = histcounts2(obj.uid,...
                obj.spikes,...
                uid_edge,...
                obj.first_event:ds:obj.last_event);
        end
        
        function n_cells_ = n_cells(obj)
            n_cells_= length(unique(obj.uid));
        end
        
        function ids_ = ids(obj)
            ids_= unique(obj.uid);
        end
        
        function duration_ = duration(obj)
            duration_ = obj.last_event - obj.first_event;
        end
        
        function disp(obj)
            obj_duration = seconds(obj.duration);
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
            
            fprintf('<%s %d units> of length %s %s \n',...
                "SpikeArray:",...
                obj.n_cells,...         
                duration_str,...
                units);
        end
        
        function n_spikes_ = n_spikes(obj)
            n_spikes_ = [];
            for ids_ = reshape(ids(obj),1,[])
                n_spikes_ = [n_spikes_;sum(obj.uid == ids_)];
            end
        end
        
        function first_event_ = first_event(obj)
            first_event_ = obj.spikes(1);
        end
        
        function last_event_ = last_event(obj)
            last_event_ = obj.spikes(end);
        end
        
        function sorted = issorted(obj)
            sorted = all(sort(obj.spikes) == sort(obj.spikes));
        end
        
        function isempty_ = isempty(obj)
            isempty_ = isempty(obj.spikes);
        end
        
        function ax = plot(obj,varargin)
            if isempty(varargin)
                ax = RasterPlot([obj.spikes,obj.uid],1);
            else
                ax = RasterPlot([obj.spikes,obj.uid],1,varargin{:});
            end
        end
    end
end