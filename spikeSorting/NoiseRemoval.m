function NoiseRemoval(basepath)
%must have already computed lfp
% come back and fancify parameters

session = bz_getSession('basepath',basepath);
nChannels = session.extracellular.nChannels;
sf = session.extracellular.srLfp;
dat_path = [basepath,filesep, session.general.name, '.dat'];
lfp = bz_GetLFP('all','basepath',basepath);


% Identify Noise intervals
mean_lfp = mean(lfp.data,2);

[bb,aa] = butter(4,[62, 120]/sf*2, 'bandpass');
hi_lfp = FilterM(bb,aa,mean_lfp);
hi_pow = abs(hilbert((hi_lfp)));
%kern = gausswin(round(.05*sf));
%hi_pow_smoo = conv(hi_pow,kern,'same');
%thresh = mean(hi_pow_smo) + 6 * std(hi_pow_smo);
%above_thresh = hi_pow_smoo>thresh;
%thresh = prctile(hi_pow,99);

thresh = mean(hi_pow) + 9 * std(hi_pow);
above_thresh = hi_pow>thresh;

clear lfp;

noiseIntsSlo = findIntervals(above_thresh);

noiseIntsFast = noiseIntsSlo*8; %adjust for 8x sampling rate in dat file [FIX THIS]

%%
% Replace in dat file
% copied from cleanPulses by manu
m = memmapfile(dat_path, 'Format','int16','Writable',true);
nPnt = round(length(m.Data)/nChannels);
intWindow = 5;

for chIdx = 1:nChannels

    for intIdx = 1:size(noiseIntsFast,1)
        noiseInds = noiseIntsFast(intIdx,1):noiseIntsFast(intIdx,2);
        frameInds = (noiseIntsFast(intIdx,1)-intWindow):(noiseIntsFast(intIdx,2)+intWindow);
        frameInds(ismember(frameInds,noiseInds)) = [];
        frameInds(frameInds>nPnt|frameInds<1) = []; 
        frameIndsFlat = sub2ind([nChannels,nPnt],chIdx*ones(size(frameInds)),frameInds);
        noiseIndsFlat = sub2ind([nChannels,nPnt],chIdx*ones(size(noiseInds)),noiseInds);
        sig = m.Data(frameIndsFlat);
        fillVals = interp1(frameInds,double(sig),noiseInds);
        m.Data(noiseIndsFlat) = int16(fillVals);
    end

end

end
