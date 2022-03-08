function synchronized = fastPETH(spiketimes,events,interval)

try
    if length(spiketimes)*length(events)<10000,
        synchronized =mainFastPETH(spiketimes,events,interval,1,length(spiketimes),length(events));
    elseif log10(length(spiketimes)*length(events))<10 ,
        synchronized = mainFastPETH(spiketimes,events,interval,10,length(spiketimes),length(events));
    else
        error('error') % in case of matrix too large
    end
catch % but also any other error
    synchronized = Sync(spiketimes,events,'durations',interval);
end
