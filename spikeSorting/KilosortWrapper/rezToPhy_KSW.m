function [spikeTimes, clusterIDs, amplitudes, templates, templateFeatures, ...
    templateFeatureInds, pcFeatures, pcFeatureInds] = rezToPhy_KSW(rez,savepath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)

% spikeTimes will be in samples, not seconds
if exist('savepath')~=1
    savepath = rez.ops.savepath;
end
fs = dir(fullfile(savepath, '*.np'));
for i = 1:length(fs)
   delete(fullfile(savepath, fs(i).name)); 
end
if exist(fullfile(savepath, '.phy'), 'dir')
    rmdir(fullfile(savepath, '.phy'), 's');
end

spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end
amplitudes = rez.st3(:,3);

Nchan = rez.ops.Nchan;

% try
%     load(rez.ops.chanMap);
% catch
%    chanMap0ind  = [0:Nchan-1]';
%    connected    = ones(Nchan, 1);
%    xcoords      = ones(Nchan, 1);
%    ycoords      = (1:Nchan)';
% end
% chanMap0 = chanMap(connected>1e-6);

connected   = rez.connected(:);
xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;

% for i = 1:length(chanMap0)
%     chanMap0(i) = chanMap0(i) - sum(chanMap0(i) > chanMap(connected<1e-6));
% end
% [~, invchanMap0] = sort(chanMap0);

templates = zeros(Nchan, nt0, rez.ops.Nfilt, 'single');
for iNN = 1:rez.ops.Nfilt
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

if ~isempty(savepath)
    
    writeNPY(spikeTimes, fullfile(savepath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1), fullfile(savepath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savepath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savepath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savepath, 'amplitudes.npy'));
    writeNPY(templates, fullfile(savepath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savepath, 'templates_ind.npy'));
    
%     Fs = rez.ops.fs;
    conn        = logical(connected);
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind(conn), fullfile(savepath, 'channel_map.npy'));
    %writeNPY(connected, fullfile(savePath, 'connected.npy'));
%     writeNPY(Fs, fullfile(savePath, 'Fs.npy'));
    writeNPY([xcoords(conn) ycoords(conn)], fullfile(savepath, 'channel_positions.npy'));
    % Added by Peter Petersen for an extra collumn in Phy with shank id
    
    writeNPY(rez.ops.kcoords, fullfile(savepath, 'channel_shanks.npy'));
    % writeNPY(rez.ops.kcoords(rez.connected(:)), fullfile(rez.ops.savepath, 'channel_shanks.npy'));
    
    writeNPY(templateFeatures, fullfile(savepath, 'template_features.npy'));
    writeNPY(templateFeatureInds'-1, fullfile(savepath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatures, fullfile(savepath, 'pc_features.npy'));
    writeNPY(pcFeatureInds'-1, fullfile(savepath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    whiteningMatrix = rez.Wrot/200;
    whiteningMatrixInv = whiteningMatrix^-1;
    writeNPY(whiteningMatrix, fullfile(savepath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savepath, 'whitening_mat_inv.npy'));
    
    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savepath, 'similar_templates.npy'));
    end
    
     %make params file
    if ~exist(fullfile(savepath,'params.py'),'file')
        fid = fopen(fullfile(savepath,'params.py'), 'w');
        
        [~, fname, ext] = fileparts(rez.ops.fbinary);
        if strcmp(savepath,rez.ops.root) 
            fprintf(fid,['dat_path = ''', fname ext, '''\n']);
        else 
            fprintf(fid,['dat_path = ''../', fname ext, '''\n']); % Peter: Added '..\' to the path to fit the custom Kilosort folder structure
        end
        fprintf(fid,['dir_path = ''''\n']); % Added by Peter Petersen
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez.ops.fs,1)
            fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
        end
        fprintf(fid,'hp_filtered = False');
        fclose(fid);
    end
end
