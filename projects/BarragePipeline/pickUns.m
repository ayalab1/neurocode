function [avgNspk, finSpkThresh, itNum, bound] = pickUns(bound, normSpkThresh, it, useMet, loadPath, Hz, ft)
basepath = pwd; basename = basenameFromBasepath(basepath);
savePath = convertStringsToChars(strcat(basepath, '\Barrage_Files\', basename, '.'));
avgNspk = 0;
itNum = 1;

if nargin < 5
    loadPath = [];
end

checking = burstCellDetection(normSpkThresh,loadPath);
while sum(checking(:,3)) < 1
    normSpkThresh = normSpkThresh - 0.1;
    checking = burstCellDetection(normSpkThresh,loadPath);
end


if it > 0
    while avgNspk <= bound
        checking = burstCellDetection(normSpkThresh,loadPath);
        load([savePath 'brstDt.cellinfo.mat']);
        
        %% We have to actually detect events, huh
        nSigma = useMet.nSigma;
        tSmooth = useMet.tSmooth;
        binSz = useMet.binSz;
        tSepMax = useMet.tSepMax; %was playing with 0.005
        mindur = useMet.mindur;
        maxdur = useMet.maxdur;
        lastmin = useMet.lastmin;
        sstd = useMet.sstd;
        estd = useMet.estd;
        EMGThresh = useMet.EMGThresh;
        neuro2 = false;
        futEVT = false;
        %With the new thresholds, we're getting it running down to UID=[],
        %which means no spikes which means no detection... need to prevent
        %this/avoid in general
        while length(spikes.UID) < 1
            bound = bound - 0.1; 
            unitsForDetection(Hz,ft);
            load([savePath 'brstDt.cellinfo.mat']);
            if bound < 1.2
                error('Bound is too low, consider manually inputting units or adjusting detection parameters');
            end
        end 
            HSE = find_HSE_b('spikes',spikes,...
                        'nSigma',nSigma,'binSz',binSz,'tSmooth',tSmooth,'tSepMax',tSepMax,...
                        'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                        'EMGThresh',EMGThresh,'sstd',sstd,'estd',estd,...
                        'neuro2',neuro2,'recordMetrics',false, 'futEVT',futEVT);        
            
        %% Get event metrics
        barSpk = getRipSpikes('basepath',pwd,'events',HSE.timestamps,'spikes',spikes,'padding',0,'saveMat',false);
        for unit = 1:length(barSpk.UnitAbs)
            for bar = 1:length(barSpk.EventAbs)
                nSpkEach(unit,bar) = length(barSpk.UnitEventAbs{unit,bar});
            end
        end

        %% Pull average number of spikes per unit per active event
        avgSpkUn = NaN(size(nSpkEach,1),1);
        for i = 1:size(nSpkEach,1)
            avgSpkUn(i) = sum(nSpkEach(i,:));
            tempCnt = length(find(nSpkEach(i,:)~=0));
            if tempCnt
                avgSpkUn(i) = avgSpkUn(i)/tempCnt;
            else
                avgSpkUn(i) = 0;
            end
        end

        %% Get our metric for success
        flag = 0;
        avgNspk = 0;
        for i = 1:length(spikes.UID)
            if avgSpkUn(i) < bound
                flag = 1;
            end
            avgNspk = avgNspk + avgSpkUn(i);
        end
        avgNspk = avgNspk/length(spikes.UID);
        if flag
            normSpkThresh = normSpkThresh + it;
            avgNspk = 0;
            itNum = itNum+1;
        end
    end
elseif it < 0
    prevThresh = 0;
    avgNspk = 1000;
    prevAvgNspk = 1000;
    keepLoop = 1;
    while (avgNspk >= bound)&&(keepLoop)
        burstCellDetection(normSpkThresh,loadPath);
        load([savePath 'brstDt.cellinfo.mat']);
        
        %% We have to actually detect events, huh
        nSigma = useMet.nSigma;
        tSmooth = useMet.tSmooth;
        binSz = useMet.binSz;
        tSepMax = useMet.tSepMax; %was playing with 0.005
        mindur = useMet.mindur;
        maxdur = useMet.maxdur;
        lastmin = useMet.lastmin;
        sstd = useMet.sstd;
        estd = useMet.estd;
        EMGThresh = useMet.EMGThresh;
        save_evts = false;
        neuro2 = false;
        futEVT = false;

        HSE = find_HSE_b('spikes',spikes,...
                    'nSigma',nSigma,'binSz',binSz,'tSmooth',tSmooth,'tSepMax',tSepMax,...
                    'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                    'EMGThresh',EMGThresh,'sstd',sstd,'estd',estd,...
                    'save_evts',save_evts,'neuro2',neuro2,'recordMetrics',false,'futEVT',futEVT);

        %% Get event metrics
        barSpk = getRipSpikes('basepath',pwd,'events',HSE.timestamps,'spikes',spikes,'padding',0,'saveMat',false);
        for unit = 1:length(barSpk.UnitAbs)
            for bar = 1:length(barSpk.EventAbs)
                nSpkEach(unit,bar) = length(barSpk.UnitEventAbs{unit,bar});
            end
        end

        %% Pull average number of spikes per unit per active event
        avgSpkUn = NaN(size(nSpkEach,1),1);
        for i = 1:size(nSpkEach,1)
            avgSpkUn(i) = sum(nSpkEach(i,:));
            tempCnt = length(find(nSpkEach(i,:)~=0));
            if tempCnt
                avgSpkUn(i) = avgSpkUn(i)/tempCnt;
            else
                avgSpkUn(i) = 0;
            end
        end

        %% Get our metric for success
        flag = 0;
        avgNspk = 0;
        for i = 1:length(spikes.UID)
            if avgSpkUn(i) < bound
                flag = 1;
            end
            avgNspk = avgNspk + avgSpkUn(i);
        end
        avgNspk = avgNspk/length(spikes.UID);
        if ~flag
            prevThresh = normSpkThresh;
            normSpkThresh = normSpkThresh + it;
            prevAvgNspk = avgNspk;
            avgNspk = 0;
            itNum = itNum+1;
        else
            normSpkThresh = prevThresh;
            avgNspk = prevAvgNspk;
            keepLoop = 0; %break
        end
    end
else
    error('You cant go anywhere! Change it from being 0');
end
finSpkThresh = normSpkThresh;

end