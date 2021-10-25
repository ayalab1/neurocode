function [spkhist spkmean spkstd ts] = spkRtHist(allspk, tSmooth, binsz)
    %% Get spike rate over time (primarily used for HSE_b)
    %quick fix for if we remove some units based on region or cell type prior to loading in
    ts = 0:binsz:allspk(end); %bins for sorting spikes
    spkhist = hist(allspk,ts); %sort spike times into those bins we created
    spkmean = mean(spkhist); %find average #spikes per bin?
    spkstd = std(spkhist);
    
    fsize = tSmooth/binsz; 
    gfilt = fspecial('gaussian',[10*fsize 1],fsize); %create our filter
    spkhist = conv(spkhist,gfilt,'same'); %filter our spikes
    spkhist = zscore(spkhist); %normalize spike times --need to think about when I should be zscoring my stuff..
    
    % addParameter(p,'tSmooth',.015,@isnumeric); % in s
    % addParameter(p,'binsz',.001,@isnumeric); % in s, originally 0.001
end