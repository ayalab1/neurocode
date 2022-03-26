function [DS1, DS2] = DetectDSpikes_v4(fname,ch_hilus,ch_molecular,res_per,varargin)

%
% 20150419 Yuta Senzai v3: consider the post-DS1 value in order to get rid of light evoked LFP deflection (for YM33)
%                          this effect seems to be only in molecular layer. This may be related to
%                          granule cell inhibition in POMC-Cre::Arch animal
%
% 20150602 Yuta Senzai v4: add m_minValue >
% res_per = default [0 Inf]
%
% INPUTS
% ch_hilus, ch_molecular: 0-based a la neuroscope
%
% default
show = 'on';

DS2_lowThreshold = 3000; % 3000*0.38uV = 1.14mV
DS1_lowThreshold = 3000; % 3000*0.38uV = 1.14mV
DS2fil_highThreshold = 3000;
DS1fil_highThreshold = 1500;
DS1wb_highThreshold = 3000;
DS2_mol_threshold = 500;

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
    end
    switch((varargin{i}))
        case 'DS2_lowThreshold'
            DS2_lowThreshold = varargin{i+1};
            if ~isdscalar(DS2_lowThreshold,'>0')
                error('Incorrect value for property ''DS2_lowThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'DS1_lowThreshold'
            DS1_lowThreshold = varargin{i+1};
            if ~isdscalar(DS1_lowThreshold,'>0')
                error('Incorrect value for property ''DS2_lowThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'DS2fil_highThreshold',
            DS2fil_highThreshold = varargin{i+1};
            if ~isdscalar(DS2fil_highThreshold,'>0'),
                error('Incorrect value for property ''DS2fil_highThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
            
        case 'DS2_mol_threshold',
            DS2_mol_threshold = varargin{i+1};
            if ~isdscalar(DS2_mol_threshold,'>0'),
                error('Incorrect value for property ''DS2_mol_threshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
            
        case 'DS1fil_highThreshold',
            DS1fil_highThreshold = varargin{i+1};
            if ~isdscalar(DS1fil_highThreshold,'>0'),
                error('Incorrect value for property ''DS1fil_highThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'DS1wb_highThreshold',
            DS1wb_highThreshold = varargin{i+1};
            if ~isdscalar(DS1wb_highThreshold,'>0'),
                error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'show',
            show = varargin{i+1};
            if ~isstring(show,'on','off'),
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
    end
end


%res_lfp_h = GetEEG(ch_hilus,'restrict',res_per);
%res_lfp_m = GetEEG(ch_molecular,'restrict',res_per);
res_lfp_h = getLFP(ch_hilus,'intervals',res_per,'noPrompts',true);
res_lfp_m = getLFP(ch_molecular,'intervals',res_per,'noPrompts',true);

fil_res_lfp_h = bz_Filter(res_lfp_h,'passband',[2 50]);
fil_res_lfp_m = bz_Filter(res_lfp_m,'passband',[2 50]);

res_lfp_h = [res_lfp_h.timestamps,double(res_lfp_h.data)];
res_lfp_m = [res_lfp_m.timestamps,double(res_lfp_m.data)];

fil_res_lfp_h = [fil_res_lfp_h.timestamps,double(fil_res_lfp_h.data)];
fil_res_lfp_m = [fil_res_lfp_m.timestamps,double(fil_res_lfp_m.data)];

% fil_res_lfp_h = FilterLFP(res_lfp_h,...
%     'nyquist',625,'passband',[2 50]);
% fil_res_lfp_m = FilterLFP(res_lfp_m,...
%     'nyquist',625,'passband',[2 50]);

time = fil_res_lfp_h(:,1);
hm_dif = fil_res_lfp_h(:,2) - fil_res_lfp_m(:,2);
%time = fil_res_lfp_h.timestamps;
%hm_dif = fil_res_lfp_h.data - fil_res_lfp_m.data;
lfp_dif = [time, hm_dif];


%% DS2 detection

DS2_thresholded = hm_dif > DS2_lowThreshold;
DS2_start = find(diff(DS2_thresholded)>0);
DS2_stop = find(diff(DS2_thresholded)<0);
% Exclude last DS if it is incomplete
if length(DS2_stop) == length(DS2_start)-1
    DS2_start = DS2_start(1:end-1);
end
% Exclude first DS if it is incomplete
if length(DS2_stop)-1 == length(DS2_start)
    DS2_stop = DS2_stop(2:end);
end
% Correct special case when both first and last DS are incomplete
if DS2_start(1) > DS2_stop(1),
    DS2_stop(1) = [];
    DS2_start(end) = [];
end

DS2_firstPass = [DS2_start,DS2_stop];

% Take out DS which are too close to file beginning- or end
DS2_firstPass = DS2_firstPass(DS2_start > 46 & DS2_stop < size(res_lfp_h,1)-45,:);

if isempty(DS2_firstPass)
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ' num2str(length(DS2_firstPass)) ' putative DS2 events.']);
end

DS2_secondPass = DS2_firstPass;


% Discard ripples with a peak power < highThreshold for
DS2 = [];
for i = 1:size(DS2_secondPass,1)
    if DS2_secondPass(i,1)>50
        m_meanValue = mean(res_lfp_m([DS2_secondPass(i,1):DS2_secondPass(i,2)],2));
        ds2criteria = mean(res_lfp_m([(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20) ],2))-DS2_mol_threshold; % 25/1250sec = 20ms, 40/1250 = 32 ms
        m_minValue =  min(res_lfp_m([(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20), (DS2_secondPass(i,2)+20):(DS2_secondPass(i,2)+45)],2));
        h_maxValue =  max(res_lfp_h([(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20), (DS2_secondPass(i,2)+20):(DS2_secondPass(i,2)+45)],2));
        %         if (df_maxValue > DS2fil_highThreshold) && (m_meanValue < ds2criteria)
        if (m_meanValue < ds2criteria && m_minValue > -6000 && h_maxValue < 6000)
            %             DS2lfp = Restrict(res_lfp_h,time(DS2_secondPass(i,:))');
            DS2lfp = res_lfp_h( DS2_secondPass(i,1): DS2_secondPass(i,2),:);
            [~,idx]=max(DS2lfp(:,2));
            DS2 = [DS2 ; time(DS2_secondPass(i,1)), DS2lfp(idx,1), time(DS2_secondPass(i,2))];
        end
    end
end

if isempty(DS2),
    disp('DS2 Peak thresholding failed.');
    return
else
    disp(['After peak thresholding: DS2 ' num2str(length(DS2)) ' events.']);
end




%% Detect DS1


DS1_thresholded = fil_res_lfp_h(:,2) > DS1_lowThreshold;
DS1_start = find(diff(DS1_thresholded)>0);
DS1_stop = find(diff(DS1_thresholded)<0);
% Exclude last DS if it is incomplete
if length(DS1_stop) == length(DS1_start)-1,
    DS1_start = DS1_start(1:end-1);
end
% Exclude first DS if it is incomplete
if length(DS1_stop)-1 == length(DS1_start),
    DS1_stop = DS1_stop(2:end);
end
% Correct special case when both first and last DS are incomplete
if DS1_start(1) > DS1_stop(1),
    DS1_stop(1) = [];
    DS1_start(end) = [];
end
DS1_firstPass = [DS1_start,DS1_stop];
if isempty(DS1_firstPass),
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ' num2str(length(DS1_firstPass)) ' putative DS1 events.']);
end

DS1_secondPass = DS1_firstPass;


% Discard DS1s with a peak power < highThreshold for
DS1 = [];
lfpsize = size(res_lfp_m,1);
for i = 1:size(DS1_secondPass,1)
    if DS1_secondPass(i,1)>80  && (DS1_secondPass(i,2)+55)<lfpsize % changed from 50 on Aug 19
        df_maxValue = max(hm_dif([DS1_secondPass(i,1):DS1_secondPass(i,2)]));
        h_maxValue = max(res_lfp_h([DS1_secondPass(i,1):DS1_secondPass(i,2)],2));
        m_meanValue = mean(res_lfp_m([DS1_secondPass(i,1):DS1_secondPass(i,2)],2));
        ds1criteria = mean(res_lfp_m([(DS1_secondPass(i,1)-80):(DS1_secondPass(i,1)-55) ],2))+500; % 25/1250sec = 20ms, 55/1250sec = 44ms
        m_maxValue = max(res_lfp_m([DS1_secondPass(i,1)-55:DS1_secondPass(i,2)+55],2));
        m_minValue = min(res_lfp_m([DS1_secondPass(i,1)-55:DS1_secondPass(i,2)+55],2));
        % make sure period does not exceed the size of res_lfp
        if size(res_lfp_h,1)>=(DS1_secondPass(i,2)+125)
            h_postValue = mean(res_lfp_h([(DS1_secondPass(i,1)+100):(DS1_secondPass(i,2)+125)],2)); % 100/1250sec = 80ms, 125/1250sec = 100ms
            ds1postcriteria = mean(res_lfp_h([(DS1_secondPass(i,1)-80):(DS1_secondPass(i,1)-55) ],2))+2000; % 25/1250sec = 20ms, 55/1250sec = 44ms
            
            
            if (df_maxValue > DS1fil_highThreshold) && (h_maxValue > DS1wb_highThreshold) && (m_meanValue > ds1criteria) ...
                    && (h_postValue < ds1postcriteria) && (m_minValue > -3000) && (m_maxValue < 3000)
                %             DS1lfp = Restrict(res_lfp_h,time(DS1_secondPass(i,:))');
                DS1lfp = res_lfp_h( DS1_secondPass(i,1): DS1_secondPass(i,2),:);
                [~,idx]=max(DS1lfp(:,2));
                DS1 = [DS1 ; time(DS1_secondPass(i,1)), DS1lfp(idx,1), time(DS1_secondPass(i,2))];
            end
        end
    end
end
if isempty(DS1)
    disp('DS1 Peak thresholding failed.');
    return
else
    disp(['After peak thresholding: DS1 ' num2str(length(DS1)) ' events.']);
end

DS1triad = DS1;
DS2triad = DS2;

% %% Synchronization
% % 'Synchronization' and 'Plotting' run very slow and are not strictly
% % necessary. To expedite, the DS2 for loop should be replaced; TH 200921
%
% % DS1
% [sync_h1,indices_h1] = Sync(res_lfp_h,DS1triad(:,2),'durations', [-0.2 0.2]);
% [sync_m1,indices_m1] = Sync(res_lfp_m,DS1triad(:,2),'durations', [-0.2 0.2]);
% nSync1 = max(indices_h1);
%
% %DS1_sync_h_population=zeros(500,nSync1);
% %DS1_sync_m_population=zeros(500,nSync1);
% DS1_sync_h_population=zeros(500,1000);
% DS1_sync_m_population=zeros(500,1000);
% DS1_goodDSidx = true(1,nSync1);
%
% for i = 1:nSync1,
%     if length(find((indices_h1==i)==1)) <500,
%         DS1_goodDSidx(i)=false;
%     elseif i<1000
%         trial1 = indices_h1==i;
%         sync_h_each = sync_h1(trial1,2);
%         sync_m_each = sync_m1(trial1,2);
%         DS1_sync_h_population(1:500,i) = sync_h_each(1:500);
%         DS1_sync_m_population(1:500,i) = sync_m_each(1:500);
%     else
%         trial1 = indices_h1==i;
%     end
% end
% DS1triad = DS1triad(DS1_goodDSidx,:);
% DS1_sync_h_population = DS1_sync_h_population(:,DS1_goodDSidx(1:1000));
% DS1_sync_m_population = DS1_sync_m_population(:,DS1_goodDSidx(1:1000));
% % DS1_sync_h_population = DS1_sync_h_population(:,DS1_goodDSidx);
% % DS1_sync_m_population = DS1_sync_m_population(:,DS1_goodDSidx);
%
% tt1 = sync_m1(trial1,1);
%
% % DS2
% [sync_h2,indices_h2] = Sync(res_lfp_h,DS2triad(:,2),'durations', [-0.2 0.2]);
% [sync_m2,indices_m2] = Sync(res_lfp_m,DS2triad(:,2),'durations', [-0.2 0.2]);
%
% % sync_h_population=[];
% % sync_m_population=[];
% nSync2 = max(indices_h2);
% % i_shift=0;
% % if length(find((indices_h2==1)==1)) < 500, nSync2 =nSync2-1; i_shift=1; DS2triad=DS2triad(2:end,1:3); end
% % if length(find((indices_h2==nSync2)==1)) < 500, nSync2 =nSync2-1; DS2triad=DS2triad(1:(end-1),1:3);  end
% % DS2_sync_h_population=zeros(500,nSync2);
% % DS2_sync_m_population=zeros(500,nSync2);
% DS2_sync_h_population=zeros(500,1000);
% DS2_sync_m_population=zeros(500,1000);
% DS2_goodDSidx = true(1,nSync2);
%
% % THIS LOOP IS EXTREMELY SLOW WITH nDS2 >> 5000; TH 200921
% for i = nSync2:-1:1%1:nSync2
%     if length(find((indices_h2==i)==1)) <500
%         DS2_goodDSidx(i)=false;
%     elseif i<1000
%         %     trial2 = indices_h2==(i+i_shift);
%         trial2 = indices_h2==i;
%         sync_h_each = sync_h2(trial2,2);
%         sync_m_each = sync_m2(trial2,2);
%         %     sync_h_population = [sync_h_population, sync_h_each(1:500)];
%         %     sync_m_population = [sync_m_population, sync_m_each(1:500)];
%         DS2_sync_h_population(1:500,i) = sync_h_each(1:500);
%         DS2_sync_m_population(1:500,i) = sync_m_each(1:500);
%     else
%         trial2 = indices_h2==i;
%     end
% end
% DS2triad = DS2triad(DS2_goodDSidx,:);
% % DS2_sync_h_population = DS2_sync_h_population(:,DS2_goodDSidx);
% % DS2_sync_m_population = DS2_sync_m_population(:,DS2_goodDSidx);
% DS2_sync_h_population = DS2_sync_h_population(:,DS2_goodDSidx(1:1000));
% DS2_sync_m_population = DS2_sync_m_population(:,DS2_goodDSidx(1:1000));
%
%
% tt2 = sync_m2(trial2,1);
%
%
% %% visualization
% %
% % if strcmp(show,'on'),
% %     figure;
% %     PlotXY(res_lfp_h,'color','r');hold on;
% %     PlotXY(res_lfp_m,'color','g');hold on;
% %     PlotXY(lfp_dif,'color','b');hold on;
% %     PlotIntervals(DS1triad(:,[1 3]), 'rectangles');hold on;
% %     PlotIntervals(DS2triad(:,[1 3]), 'rectangles');hold on;
% % end
%
% % DS1
% if strcmp(show,'on'),
%     fig1=figure;
%     plot(tt1(1:500),DS1_sync_h_population','r');hold on;
%     plot(tt1(1:500),DS1_sync_m_population','g');hold on;
%     %SaveFigPng(fig1,[fname '--DS2-LFP-individual-hilusCh' num2str(ch_hilus) '-moleCh' num2str(ch_molecular)]);
% %     fig2=figure;
% %     shadedErrorBar(tt1(1:500), DS1_sync_h_population', {@mean,@std},'r',1);hold on;
% %     shadedErrorBar(tt1(1:500), DS1_sync_m_population', {@mean,@std},'g',1);hold on;
% %     SaveFigPng(fig2,[fname '--DS1-LFP-populstion-hilusCh' num2str(ch_hilus) '-moleCh' num2str(ch_molecular)]);
% end
%
% % DS2
% if strcmp(show,'on'),
%     fig3=figure;
%     plot(tt2(1:500),DS2_sync_h_population','r');hold on;
%     plot(tt2(1:500),DS2_sync_m_population','g');hold on;
%     %SaveFigPng(fig3,[fname '--DS2-LFP-individual-hilusCh' num2str(ch_hilus) '-moleCh' num2str(ch_molecular)]);
% %     fig4=figure;
% %     shadedErrorBar(tt2(1:500), DS2_sync_h_population', {@mean,@std},'r',1);hold on;
% %     shadedErrorBar(tt2(1:500), DS2_sync_m_population', {@mean,@std},'g',1);hold on;
% %     SaveFigPng(fig4,[fname '--DS2-LFP-population-hilusCh' num2str(ch_hilus) '-moleCh' num2str(ch_molecular)]);
% end

%% save events
% .evt (FMA standard)
n = size(DS1triad,1);
d1 = DS1triad(:,1:3)';
events1.time = d1(:);
for i = 1:3:3*n,
    events1.description{i,1} = ['DentateSpike1 start ' num2str(ch_hilus)];
    events1.description{i+1,1} = ['DentateSpike1 peak ' num2str(ch_hilus) ];
    events1.description{i+2,1} = ['DentateSpike1 stop ' num2str(ch_hilus) ];
end

n = size(DS2triad,1);
d2 = DS2triad(:,1:3)';
events2.time = d2(:);
for i = 1:3:3*n,
    events2.description{i,1} = ['DentateSpike2 start ' num2str(ch_hilus)];
    events2.description{i+1,1} = ['DentateSpike2 peak ' num2str(ch_hilus)];
    events2.description{i+2,1} = ['DentateSpike2 stop ' num2str(ch_hilus)];
end

% save([fname '.DSsynch.ch' num2str(ch_hilus) '.mat'],'DS1_sync_h_population','DS1_sync_m_population','DS2_sync_h_population','DS2_sync_m_population','tt1','tt2');

% SaveEvents([fname '-hilus_ch' num2str(ch_hilus) '-ml_ch' num2str(ch_molecular) '.DS1.evt'],events1);
% SaveEvents([fname '-hilus_ch' num2str(ch_hilus) '-ml_ch' num2str(ch_molecular) '.DS2.evt'],events2);

% .event.mat (buzcode & cell_explorer standard)
clear DS1 DS2

DS1.timestamps = DS1triad(:,[1 3]);
DS1.peaks = DS1triad(:,2);
DS1peakIDCs = round(DS1triad(:,2)*1250);
DS1.amplitudes = (res_lfp_h(DS1peakIDCs,2)-res_lfp_m(DS1peakIDCs,2))*.00038; % .00038 = conversion to mV
% amplitude: amplitude of each event (Px1).
% amplitudeUnits: specify the units of the amplitude vector.
DS1.amplitudeUnits = 'mV';
% eventID: numeric ID for classifying various event types (Px1).
DS1.eventID = ones(size(DS1triad,1),1);
% eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
DS1.eventIDlabels = repmat({'DS1'},size(DS1triad,1),1);
% eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
DS1.eventIDbinary = false(size(DS1triad,1),1);
% center: center time-point of event (in seconds; calculated from timestamps; Px1).
% duration: duration of event (in seconds; calculated from timestamps; Px1).
DS1.duration = DS1triad(:,2)-DS1triad(:,1);
DS1.center = DS1triad(:,1)+DS1.duration;

DS1.detectorinfo.detectorname = 'DetectDSpikes_v4';
DS1.detectorinfo.detectionparms = [];
DS1.detectorinfo.detectionintervals = [0 Inf];
DS1.detectorinfo.detectiondate = datetime('today');
DS1.detectorinfo.ml_channel = ch_molecular;
DS1.detectorinfo.h_channel = ch_hilus;

DS2.timestamps = DS2triad(:,[1 3]);
DS2.peaks = DS2triad(:,2);
DS2peakIDCs = round(DS2triad(:,2)*1250);
DS2.amplitudes = (res_lfp_h(DS2peakIDCs,2)-res_lfp_m(DS2peakIDCs,2))*.00038; % .00038 = conversion to mV
% amplitude: amplitude of each event (Px1).
% amplitudeUnits: specify the units of the amplitude vector.
DS2.amplitudeUnits = 'mV';
% eventID: numeric ID for classifying various event types (Px1).
DS2.eventID = ones(size(DS2triad,1),1);
% eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
DS2.eventIDlabels = repmat({'DS2'},size(DS2triad,1),1);
% eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
DS2.eventIDbinary = false(size(DS2triad,1),1);
% center: center time-point of event (in seconds; calculated from timestamps; Px1).
% duration: duration of event (in seconds; calculated from timestamps; Px1).
DS2.duration = DS2triad(:,2)-DS2triad(:,1);
DS2.center = DS2triad(:,1)+DS2.duration;

DS2.detectorinfo.detectorname = 'DetectDSpikes_v4';
DS2.detectorinfo.detectionparms = [];
DS2.detectorinfo.detectionintervals = [0 Inf];
DS2.detectorinfo.detectiondate = datetime('today');
DS2.detectorinfo.ml_channel = ch_molecular;
DS2.detectorinfo.h_channel = ch_hilus;

save([fname '.DS1.events.mat'],'DS1');
save([fname '.DS2.events.mat'],'DS2');

end