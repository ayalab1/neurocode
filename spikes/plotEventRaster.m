function  plotEventRaster(event,varargin)
%
%   Plot spike raster for invidual events (such as ripples), sorting cells
%   by firing order and color code them according to diverse features
%   (region, cell type, etc.)
%
%   Inputs:
%   event   = [start stop] in seconds for one or multiple events
%   spikes  = structure of spike time and UID info
%   lfpChan = channel to plot lfp (base 1)
%   loadDat = load lfp trace from .dat instead of .lfp. Default = false
%   tag     = feature to color code raste. Now supporting: pyrInt, 
%               brainRegion, deepSup, REMshift
%   tag2    = feature for shape of the raster. Now supporting: cellType,
%               ripMod(**BUT ONLY IN CONJUNCTION WITH BRAIN REGION AS TAG**) 
%   saveFig = default false

%   Antonio FR, 10/21. Lindsay Karaba, 11/21

%% inputs
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'lfpChan',[],@isnumeric);
addParameter(p,'tag','brainRegion',@isstr);
addParameter(p,'tag2','cellType',@isstr);
addParameter(p,'savePath',pwd,@isstr);
addParameter(p,'evtNum',[],@isnumeric);
addParameter(p,'loadDat',false,@islogical);
addParameter(p,'regions',[],@isstring);
addParameter(p,'SleepState',[],@isstruct);
addParameter(p,'saveFig',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
spikes = p.Results.spikes;
lfpChan = p.Results.lfpChan;
tag = p.Results.tag;
tag2 = p.Results.tag2;
savePath = p.Results.savePath;
evtNum = p.Results.evtNum;
loadDat = p.Results.loadDat;
regions = p.Results.regions;
SleepState = p.Results.SleepState;
saveFig = p.Results.saveFig;

basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
load(fullfile(basepath,[basename '.session.mat']));
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']));

if isempty(spikes)
    load(fullfile(basepath,[basename '.spikes.cellinfo.mat']));
end

if strcmp(tag2,'ripMod')&&~strcmp(tag,'brainRegion')
    error('ripMod is not implemented with anything other than brainRegion color coding');
end

sr = session.extracellular.sr;

pad = 0.05; % padding time
for e = 1:size(event,1)
    event(e,1) = event(e,1)-pad;
    event(e,2) = event(e,2)+pad;
end

if ~isempty(SleepState)
	for e = 1:size(event,1)
		stateName{e} = getCurState(SleepState, event(e,:));
	end
end

%% plot lfp
if ~isempty(lfpChan)
    lfp = [];
    if ~loadDat
        lfp = getLFP(lfpChan,'intervals',event,'basepath',basepath);
        % add option to filter LFP
        figure('Position', get(0, 'Screensize'));  
        for e = 1:size(event,1)
            subplot(2,size(event,1),e);
            plot(lfp(e).timestamps,lfp(e).data(:,1),'k');hold on;
            if ~isempty(SleepState)
				title(stateName{e});
			end
            xlim([lfp(e).timestamps(1) lfp(e).timestamps(end)]);
        end
    elseif loadDat
        figure('Position',[800 400 900 500]);  
        for e = 1:size(event,1)
            lfpdat = LoadBinary([basename '.dat'],'frequency',sr,'nChannels',session.extracellular.nChannels,...
                'channels',lfpChan,'start',event(e,1),'duration',event(e,2)-event(e,1));  
            lfp(e).data = lfpdat; 
            t = event(e,1):(1/sr):event(e,2);   lfp(e).timestamps = t(1:size(lfp(e).data(:,1),1))'; %sometimes the time array doesn't match?
            subplot(2,size(event,1),e);
            plot(lfp(e).timestamps,lfp(e).data(:,1),'k');hold on;
            if ~isempty(SleepState)
				title(stateName{e});
			end
            xlim([lfp(e).timestamps(1) lfp(e).timestamps(end)]);
        end
    end
end

%% plot spike raster: pyr int
for e = 1:size(event,1)
    t = lfp(e).timestamps; % add alternative for no lfp
    
    rasterT = zeros(length(spikes.times),length(t));
    for i = 1:length(spikes.times) % colect spk in rip
        temp{i} = Restrict(spikes.times{i},[t(1) t(end)]);
    end
    for i = 1:length(spikes.times)
        if ~isempty(temp{i})
            for j = 1:length(temp{i})
                [a,b] = min((abs(t-temp{i}(j)))); % b = indice of spike in t vector
                rasterT(i,b) = temp{i}(j); % matrix con ts
                if j == 1
                    firstSpk(i,1) = b;
                end
                clear a b
            end
        else
            firstSpk(i,1) = 0;
        end
    end
    [a,b]= sort(firstSpk);
    rasterO = rasterT(b,:);
    rasterIDo = spikes.UID(b)';
    
    count = 0;
    for i = 1:size(rasterO,1)
        if sum(rasterO(i,:)) > 0
            count = count +1;
            raster(count,:) = rasterO(i,:);
            rasterID(count,:) = rasterIDo(i,:);
        end
    end
    subplot(2,size(event,1),e+size(event,1));
    hold on;
    yticks([1:length(rasterID)]);
    yticklabels(num2str(rasterID(:)));
    if strcmp('cellType',tag2)
        for i = 1:size(raster,1)
            for j = 1:size(raster,2)
                if raster(i,j) > 0 && strcmp('Narrow Interneuron',cell_metrics.putativeCellType{rasterID(i)})
                    scatter(raster(i,j),i,'.k');hold on;
                    clear y;
                elseif raster(i,j) > 0  strcmp('Pyramidal Cell',cell_metrics.putativeCellType{rasterID(i)});
                    scatter(raster(i,j),i,'vk','filled');hold on;
                    clear y;
                end
            end
        end
        xlim([t(1) t(end)]);ylim([0 size(raster,1)+1]);
        clear rasterT temp firstSpk raster a b
    elseif strcmp('ripMod',tag2)
        for i = 1:size(raster,1)
            for j = 1:size(raster,2)
                if raster(i,j) > 0 && ismember(rasterID(i),cell_metrics.tags.N)
                    scatter(raster(i,j),i,'.k');hold on;
                    clear y;
                elseif raster(i,j) > 0 && ismember(rasterID(i),cell_metrics.tags.P)
                    scatter(raster(i,j),i,'vk','filled');hold on;
                    clear y;
                end
            end
        end
        xlim([t(1) t(end)]);ylim([0 size(raster,1)+1]);
        clear rasterT temp firstSpk raster a b
    else
        error('cellType and ripMod are the only two secondary tags currently implemented');
    end
    
end

%% color spikes by tag
switch(tag)
    case 'pyrInt'
        
    case 'brainRegion'
        if isempty(regions)
			regions = unique(cell_metrics.brainRegion);
			colors = distinguishable_colors(numel(regions));
		else
			colors = ['m';'g';'b';'r';'c';'k'];
			if length(colors) <= length(regions)
				colors = colors(1:length(regions));
			else
				warning('Too many regions, defaulting to distinguishable colors');
				regions = unique(cell_metrics.brainRegion);
				colors = distinguishable_colors(numel(regions));
			end
		end
        for e = 1:size(event,1)
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            t = lfp(e).timestamps; % add alternative for no lfp
            
            rasterT = zeros(length(spikes.times),length(t));
            for i = 1:length(spikes.times) % colect spk in rip
                temp{i} = Restrict(spikes.times{i},[t(1) t(end)]);
            end
            for i = 1:length(spikes.times)
                if ~isempty(temp{i})
                    for j = 1:length(temp{i})
                        [a,b] = min((abs(t-temp{i}(j)))); % b = indice of spike in t vector
                        rasterT(i,b) = temp{i}(j); % matrix con ts
                        if j == 1
                            firstSpk(i,1) = b;
                        end
                        clear a b
                    end
                else
                    firstSpk(i,1) = 0;
                end
            end
            [a,b]= sort(firstSpk);
            rasterO = rasterT(b,:);
            rasterIDo = spikes.UID(b)';clear a b;
            
            count = 0;
            for i = 1:size(rasterO,1)
                if sum(rasterO(i,:)) > 0
                    count = count +1;
                    raster(count,:) = rasterO(i,:);
                    rasterID(count,:) = rasterIDo(i,:);
                end
            end
            
            subplot(2,size(event,1),e+size(event,1));
            for i = 1:size(raster,1)
                for j = 1:size(raster,2)
                    if strcmp('cellType',tag2)
                        if raster(i,j) > 0 && strcmp('Narrow Interneuron',cell_metrics.putativeCellType{rasterID(i)})
                            br = find(strcmp(regions,cell_metrics.brainRegion(rasterID(i))));
                            if isempty(br); br = 1; end
                            scatter(raster(i,j),i,10,colors(br,:),'o','filled');hold on;
                            clear y;
                        elseif raster(i,j) > 0 && strcmp('Pyramidal Cell',cell_metrics.putativeCellType{rasterID(i)})
                            br = find(strcmp(regions,cell_metrics.brainRegion(rasterID(i))));
                            if isempty(br); br = 1; end
                            scatter(raster(i,j),i,30,colors(br,:),'v','filled');hold on;
                            clear y;
                        end
                    elseif strcmp('ripMod',tag2)
                        if raster(i,j) > 0 && ismember(rasterID(i),cell_metrics.tags.N)
                            br = find(strcmp(regions,cell_metrics.brainRegion(rasterID(i))));
                            if isempty(br); br = 1; end
                            scatter(raster(i,j),i,10,colors(br,:),'o','filled');hold on;
                            clear y;
                        elseif (raster(i,j) > 0) && (ismember(rasterID(i),cell_metrics.tags.P))
                            br = find(strcmp(regions,cell_metrics.brainRegion(rasterID(i))));
                            if isempty(br); br = 1; end
                            scatter(raster(i,j),i,30,colors(br,:),'v','filled');hold on;
                            clear y;
                        end
                    else
                        error('cellType and ripMod are the only two secondary tags currently implemented');
                    end
                end
            end
            xlim([t(1) t(end)]);ylim([0 size(raster,1)+1]);
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            
        end
        
    case 'deepSup'
        for e = 1:size(event,1)
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            t = lfp(e).timestamps; % add alternative for no lfp
            
            rasterT = zeros(length(spikes.times),length(t));
            for i = 1:length(spikes.times) % colect spk in rip
                temp{i} = Restrict(spikes.times{i},[t(1) t(end)]);
            end
            for i = 1:length(spikes.times)
                if ~isempty(temp{i})
                    for j = 1:length(temp{i})
                        [a,b] = min((abs(t-temp{i}(j)))); % b = indice of spike in t vector
                        rasterT(i,b) = temp{i}(j); % matrix con ts
                        if j == 1
                            firstSpk(i,1) = b;
                        end
                        clear a b
                    end
                else
                    firstSpk(i,1) = 0;
                end
            end
            [a,b]= sort(firstSpk);
            rasterO = rasterT(b,:);
            rasterIDo = spikes.UID(b)';clear a b;
            
            count = 0;
            for i = 1:size(rasterO,1)
                if sum(rasterO(i,:)) > 0
                    count = count +1;
                    raster(count,:) = rasterO(i,:);
                    rasterID(count,:) = rasterIDo(i,:);
                end
            end
            
            subplot(2,size(event,1),e+size(event,1));
            for i = 1:size(raster,1)
                for j = 1:size(raster,2)
                    if raster(i,j) > 0 && strcmp('Narrow Interneuron',cell_metrics.putativeCellType{rasterID(i)})
                        scatter(raster(i,j),i,10,'k','o','filled');hold on;clear y;
                    elseif raster(i,j) > 0  strcmp('Pyramidal Cell',cell_metrics.putativeCellType{rasterID(i)});
                        if strcmp('Deep',cell_metrics.deepSuperficial{rasterID(i)})
                            scatter(raster(i,j),i,30,'b','v','filled');hold on;clear y;
                        elseif strcmp('Superficial',cell_metrics.deepSuperficial{rasterID(i)})
                            scatter(raster(i,j),i,30,'r','v','filled');hold on;clear y;
                        end
                    end
                end
            end
            xlim([t(1) t(end)]);ylim([0 size(raster,1)+1]);
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            
        end
        
    case 'REMshift'
        load([basename '.theta_rem_shift.mat']);
        for e = 1:size(event,1)
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            t = lfp(e).timestamps; % add alternative for no lfp
            
            rasterT = zeros(length(spikes.times),length(t));
            for i = 1:length(spikes.times) % colect spk in rip
                temp{i} = Restrict(spikes.times{i},[t(1) t(end)]);
            end
            for i = 1:length(spikes.times)
                if ~isempty(temp{i})
                    for j = 1:length(temp{i})
                        [a,b] = min((abs(t-temp{i}(j)))); % b = indice of spike in t vector
                        rasterT(i,b) = temp{i}(j); % matrix con ts
                        if j == 1
                            firstSpk(i,1) = b;
                        end
                        clear a b
                    end
                else
                    firstSpk(i,1) = 0;
                end
            end
            [a,b]= sort(firstSpk);
            rasterO = rasterT(b,:);
            rasterIDo = spikes.UID(b)';clear a b;
            
            count = 0;
            for i = 1:size(rasterO,1)
                if sum(rasterO(i,:)) > 0
                    count = count +1;
                    raster(count,:) = rasterO(i,:);
                    rasterID(count,:) = rasterIDo(i,:);
                end
            end
            
            subplot(2,size(event,1),e+size(event,1));
            for i = 1:size(raster,1)
                for j = 1:size(raster,2)
                    if raster(i,j) > 0 && strcmp('Narrow Interneuron',cell_metrics.putativeCellType{rasterID(i)})
                        scatter(raster(i,j),i,10,'k','o','filled');hold on;clear y;
                    elseif raster(i,j) > 0  strcmp('Pyramidal Cell',cell_metrics.putativeCellType{rasterID(i)});
                        if rem_shift_data.rem_shift(rasterID(i)) == 1
                            scatter(raster(i,j),i,30,'b','v','filled');hold on;clear y;
                        elseif rem_shift_data.non_rem_shift(rasterID(i)) == 1
                            scatter(raster(i,j),i,30,'r','v','filled');hold on;clear y;
                        else
                            scatter(raster(i,j),i,30,'k','v','filled');hold on;clear y;
                        end
                    end
                end
            end
            xlim([t(1) t(end)]);ylim([0 size(raster,1)+1]);
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            
        end
        
    case 'ripParticip'
        for e = 1:size(event,1)
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
            t = lfp(e).timestamps; % add alternative for no lfp
            
            rasterT = zeros(length(spikes.times),length(t));
            for i = 1:length(spikes.times) % colect spk in rip
                temp{i} = Restrict(spikes.times{i},[t(1) t(end)]);
            end
            for i = 1:length(spikes.times)
                if ~isempty(temp{i})
                    for j = 1:length(temp{i})
                        [a,b] = min((abs(t-temp{i}(j)))); % b = indice of spike in t vector
                        rasterT(i,b) = temp{i}(j); % matrix con ts
                        if j == 1
                            firstSpk(i,1) = b;
                        end
                        clear a b
                    end
                else
                    firstSpk(i,1) = 0;
                end
            end
            [a,b]= sort(firstSpk);
            rasterO = rasterT(b,:);
            rasterIDo = spikes.UID(b)';clear a b;
            
            count = 0;
            for i = 1:size(rasterO,1)
                if sum(rasterO(i,:)) > 0
                    count = count +1;
                    raster(count,:) = rasterO(i,:);
                    rasterID(count,:) = rasterIDo(i,:);
                end
            end
            
            subplot(2,size(event,1),e+size(event,1));
            for i = 1:size(raster,1)
                for j = 1:size(raster,2)
                    if raster(i,j) > 0 && strcmp('Narrow Interneuron',cell_metrics.putativeCellType{rasterID(i)})
                        scatter(raster(i,j),i,10,'k','o','filled');hold on;clear y;
                    elseif raster(i,j) > 0  strcmp('Pyramidal Cell',cell_metrics.putativeCellType{rasterID(i)});
                        if cell_metrics.ripple_particip(rasterID(i)) > 0.4
                            scatter(raster(i,j),i,30,'b','v','filled');hold on;clear y;
                        elseif cell_metrics.ripple_particip(rasterID(i)) < 0.1
                            scatter(raster(i,j),i,30,'r','v','filled');hold on;clear y;
                        else
                            scatter(raster(i,j),i,30,'k','v','filled');hold on;clear y;
                        end
                    end
                end
            end
            xlim([t(1) t(end)]);ylim([0 size(raster,1)+1]);
            clear rasterT temp firstSpk raster rasterO rasterID rasterIDo
        end
        
    otherwise
        error(['Unknown property ']);
end

% %% add saving
if saveFig
    oldPath = cd(savePath);
    rastFiles = dir(['*' animName '.' basename '.rastEvt*.png']);
    cd(oldPath);
    if isempty(rastFiles)
        fileN = 1;
    else
        % Set file index to next available value
        pat = ['.rastEvt[0-9].'];
        fileN = 0;
        for ii = 1:length(rastFiles)
            token  = regexp(rastFiles(ii).name,pat);
            val    = str2double(rastFiles(ii).name(token+9:token+10));
            fileN  = max([fileN val]);
        end
        fileN = fileN + 1;
    end
    if fileN < 10
        useFileN = strcat('0',num2str(fileN));
    else
        useFileN = num2str(fileN);
    end
    saveas(gcf,[savePath '\' animName '.' basename '.rastEvt' useFileN '.png']);
end

end




