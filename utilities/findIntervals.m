function out = findIntervals(booString)
%findIntervals - find start and end indices from boolstring
booString = reshape(booString,[1,length(booString)]);

starts = find(diff([false booString])>0);
ends = find(diff([booString false])<0);

out = [starts', ends'];

end