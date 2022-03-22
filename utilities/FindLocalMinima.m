function troughs = FindLocalMinima(signal)

if isdmatrix(signal,'@2'),
    iTroughs = FindLocalMinima(signal(:,2));
    troughs = signal(iTroughs,1);
    return
end

d = [nan;diff(signal(:))>0];
troughs = strfind(d',[0 1])';
