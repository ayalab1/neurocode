function stateName = getCurState(SleepState, event)
%getCurState - Get current state name.
% comp = ['WAKEstate', 'NREMstate', 'REMstate', 'THETA', 'nonTHETA', 'THETAtask'];
fields = fieldnames(SleepState.ints);
state = cell(2,1);
for t = 1:size(event,2)
	for f = 1:length(fields)
		curState = fields{f};
        switch(curState)
            case 'WAKEstate'
                compTs = SleepState.ints.WAKEstate;
            case 'NREMstate'
                compTs = SleepState.ints.NREMstate;
            case 'REMstate'
                compTs = SleepState.ints.REMstate;
            case 'THETA'
                compTs = SleepState.ints.THETA;
            case 'nonTHETA'
                compTs = SleepState.ints.nonTHETA;
            case 'THETAtask'
                compTs = SleepState.ints.THETAtask;
            otherwise 
                error('Unknown sleep state, please update case conditions');
        end    
        checkStart = [];
        checkEnd = [];
        checkStart = compTs(:,1)<=event(t);
        checkStart = find(checkStart == 1);
        checkEnd = compTs(:,2)>=event(t);
        checkEnd = find(checkEnd == 1);
        [~,~,checkState] = intersect(checkStart, checkEnd);
        if ~isempty(checkState)
            if isempty(state{t})
                state{t} = convertCharsToStrings(curState);
            else
                state{t} = strcat(state{t},', ',convertCharsToStrings(curState));
            end
        end
    end
    if isempty(state{t})
        state{t} = "Unknown";
    end
end

if strcmp(state{1},state{2})
    stateName = state{1};
else
    stateName = strcat(state{1},'->',state{2});
end
				
end