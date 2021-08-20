function out = findIntervals(booString)

booString = reshape(booString,[1,length(booString)]);

starts = find(diff([false booString])>0);
ends = find(diff([booString false])<0);

out = [starts', ends'];

end