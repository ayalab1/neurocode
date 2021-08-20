function PhyAutoClustering_km(clusteringpath,varargin)
% AutoClustering automtically cleans Kilosort output in phy format defined
% by a clusteringpath.
%
% INPUT:
%     clusteringpath: char
%
% optional:
%     AutoClustering(clusteringpath,elec,dim)
%     where dim is the number of channels in electro group (if not
%     defined, will read the first line of the fet file
%
% AutoClustering is meant to clean the output of KlustaKwik. The first
% thing it does is to separate electrical artifacts and MUA from putative
% isolated units. To do so, it sorts out units which have no clear
% refractory period (based on Hill, Mehta and Kleinfeld, J Neurosci.,
% 2012). Threshold can be set in the parameter section of this file
% ("Rogue spike threshold"). Then, it separates electrical
% artifats from MUA based on the assumption that electrical artifacts are
% highly correlated on the different channels: the average waveform of at
% least one channel has to be different from the across-channel average
% waveform by a certrain amount of total variance (can be set in the
% parameter section, "Deviation from average spike threshold")
%
%
% Once the program has determined which of the clusters are putative
% isolated units, it tries to merge them based on waveform similarity
% (mahalanobis distance) and quality of the refractory period in the new
% merged cluster (or "Inter Common Spike Interval" from MS Fee et al. J
% Neurosci. Meth., 1996)
%
% Original script by Adrien Peyrache, 2012.
% Many modifications for Phy processing pipeline by 
% Yuta Senzai and Peter Petersen


% if ~isempty(varargin)
%     dim = varargin{1};
%     dim = dim(:);
%     if any(double(int16(dim))~=dim)
%         error('Number of dimensions must be an integer')
%     end
%     
%     if size(dim,1) ~= length(elec) && length(dim) ~=1
%         error('Number of dimensions must be a vector of the same length as electrode vector or a single value')
%     end
%     if length(dim) == 1
%         dim = dim*ones(length(elec),1);
%     end
% else
%     dim = zeros(length(elec),1);
% end

% Rogue spike threshold (for MUA); value between 0 an 1
%rogThres = 0.25;
rogThres = 0.33;

% Relative deviation (from 0 to 1) from average spike threshold (for electrical artifacts)
%devThres = 0.25;
% =1000 => bypass it
% devThres = 1000;
rThres = 0.7;
mprThres = 2;

% Artifact removal threshold
amplitude_thr = 50;
mahal_thr = 18;

% Load spike timing
cd(clusteringpath)
dirname = ['PhyAutoClustering_', datestr(clock,'yyyy-mm-dd_HHMMSS')];
mkdir(dirname)
copyfile(fullfile(clusteringpath, 'spike_clusters.npy'), fullfile(clusteringpath, dirname,'spike_clusters.npy'))
if exist(fullfile(clusteringpath, 'cluster_group.tsv'))
    copyfile(fullfile(clusteringpath, 'cluster_group.tsv'), fullfile(clusteringpath, dirname, 'cluster_group.tsv'))
end

clu = readNPY(fullfile(clusteringpath, 'spike_clusters.npy'));
clu = double(clu);
cids = unique(clu);

wav_all_orig = readNPY(fullfile(clusteringpath,'templates.npy'));
wav_all_orig2 = permute(wav_all_orig,[2,3,1]);
ch_indx = [];
for i = 1:size(wav_all_orig2,3)
    [~,ch_indx(i)] = max(max(wav_all_orig2(:,:,i))-min(wav_all_orig2(:,:,i)));
end
channel_shanks = readNPY(fullfile(clusteringpath, 'channel_shanks.npy'));

ch_indx2 = {};
for j = unique(channel_shanks)
    ch_indx2{j} = find(channel_shanks == j);
end

% Removing spikes with large artifacts
disp('Removing spikes with large artifacts')
spike_amplitudes = readNPY(fullfile(clusteringpath, 'amplitudes.npy'));

spike_amplitudes = nanconv(spike_amplitudes',ones(1,20),'edge'); % what do you get with the convolution????
indx = find(spike_amplitudes>amplitude_thr);
disp([num2str(length(indx)),' artifacts detected'])

clu2 = clu;
artifact_clusters = [];
for j = unique(channel_shanks)
    indx22 = find(ismember(ch_indx(clu(indx)+1),ch_indx2{j}));
    clu2(indx(indx22)) = max(clu)+j;
    artifact_clusters = [artifact_clusters,max(clu)+j];
end
clu = clu2;
spike_PCAs = double(readNPY(fullfile(clusteringpath, 'pc_features.npy')));

% Mahal artifact removal
disp('Removing outliers by Mahalanobis theshold...')
spike_clusters = clu;
mahal_outlier_clusters = [];
spikes_removed = 0;
for i = 1:length(cids)
    cluster_id = cids(i);
    indexes = find(spike_clusters==cluster_id);
    if length(indexes)>100
        indexes1 = spike_PCAs(indexes,:,:);
        indexes2 = reshape(indexes1,[size(indexes1,1),size(indexes1,2)*size(indexes1,3)]);
        test2 = mahal(indexes2,indexes2);
        test3 = find(test2>mahal_thr^2);
        mahal_outlier_clusters = [mahal_outlier_clusters,max(spike_clusters)+1];
        spikes_removed = spikes_removed+length(test3);
        spike_clusters(indexes(test3)) = mahal_outlier_clusters(end);
    end
end
clu = spike_clusters;
disp([num2str(length(mahal_outlier_clusters)),' units cleaned by Mahal outlier detection. Spikes removed: ',num2str(spikes_removed)])
writeNPY(uint32(clu), fullfile(clusteringpath, 'spike_clusters.npy'));

% % Loading rez.mat for sampling rate
% disp('Loading rez.mat')
% load(fullfile(clusteringpath,'rez.mat'))
% sr = rez.ops.fs;
% 
% res_int = readNPY(fullfile(clusteringpath,'spike_times.npy'));
% res = double(res_int)/sr;
% 
% wav_all = wav_all_orig2;
% 
% disp('Classifying noise/mua')
% meanR = [];
% fractRogue = [];
% maxPwRatio = [];
% for ii=1:length(cids)
%     spktime = res(clu==cids(ii));
%     if ~isempty(spktime) && ~any(any(isnan(squeeze(wav_all(13:end,:,ii)))))
%         %     dim = channel_shanks(ch_indx(cids(ii)+1));
%         wav = squeeze(wav_all(13:end,:,ii));
%         wav = wav(:,find(any(wav)));
%         dim = size(wav,2);
%         
%         [R,~] = corrcoef(wav);
%         meanR_cur = (sum(sum(R)) - dim) /(dim*(dim-1));
%         meanR = [meanR; meanR_cur];
%         
%         maxPwRatio_cur = max(abs(wav(11,:)))/mean(abs(wav(11,:)));
%         maxPwRatio = [maxPwRatio; maxPwRatio_cur];
%         
%         [ccgR,t] = CCG(spktime,ones(size(spktime)),'binsize',.0005,'duration',.06);
%         indx3 = find(t > -0.0015 & t < 0.0015);
%         spkRef  = mean(ccgR(indx3)); % refractory period: -1.5ms to 1.5ms
%         spkMean = mean(ccgR(round(indx3(1)/2):indx3(1)-1));
%         %         l = FractionRogueSpk(spktime,tR,tC);
%         l = spkRef/spkMean;
%         fractRogue = [fractRogue;l];
%     else
%         maxPwRatio = [maxPwRatio; 0];
%         meanR = [meanR; 0];
%         fractRogue = [fractRogue;0];
%     end
% end
% % Here we compute # of spike per cell. Some code for the errormatrix fails
% % when the cluster is defined by only a few samples. We'll put a
% % threshopld a bit later on the total # of spikes.
% h = hist(clu,unique(clu));
% h = h(:);
% h = h(1:length(meanR));
% % Definition of cluster 0 (noiseIx) and cluster 1 (muaIx)
% % Outliers of total spike power (putative electrical artifacts) not imlemented yet
% noiseIx = find((meanR >= rThres & maxPwRatio < mprThres)|h<100);
% muaIx = find(fractRogue>rogThres & ~(meanR >= rThres & maxPwRatio < mprThres) & h>=100);
% goodIx = find(fractRogue<=rogThres & ~(meanR >= rThres & maxPwRatio < mprThres) & h>=100); % 100 or samlenum
% 
% % Saving clusters to cluster_group.tsv
% fid = fopen(fullfile(clusteringpath,'cluster_group.tsv'),'w');
% fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
% for ii=1:length(cids)
%     if any(clu==cids(ii))
%         if any(goodIx==ii)
%             fwrite(fid, sprintf('%d\t%s\r\n', cids(ii), 'good'));
%         elseif any(muaIx==ii)
%             fwrite(fid, sprintf('%d\t%s\r\n', cids(ii), 'mua'));
%         elseif any(noiseIx==ii)
%             fwrite(fid, sprintf('%d\t%s\r\n', cids(ii), 'noise'));
%         end
%     end
% end
% for jj = 1:length(artifact_clusters)
%     fwrite(fid, sprintf('%d\t%s\r\n', artifact_clusters(jj), 'artifacts'));
% end
% mahal_outlier_clusters = unique(mahal_outlier_clusters);
% for jj = 1:length(mahal_outlier_clusters)
%     fwrite(fid, sprintf('%d\t%s\r\n', mahal_outlier_clusters(jj), 'mua'));
% end
% fclose(fid);
% 
% % save(fullfile(clusteringpath,'autoclusta_params.mat'),'meanR','maxPwRatio','fractRogue','noiseIx','muaIx');
% disp('AutoClustering complete.')
end
