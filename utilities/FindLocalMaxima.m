function peaks = FindLocalMaxima(signal)

if isdmatrix(signal,'@2'),
    iPeaks = FindLocalMaxima(signal(:,2));
    peaks = signal(iPeaks,1);
    return
end

d = [nan;diff(signal(:))>0];
peaks = strfind(d',[1 0])';
