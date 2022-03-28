function [DS1, DS2] = DetectDSpikes_v4(ch_hilus,ch_molecular,varargin)

%
% 20150419 Yuta Senzai v3: consider the post-DS1 value in order to get rid of light evoked LFP deflection (for YM33)
%                          this effect seems to be only in molecular layer. This may be related to
%                          granule cell inhibition in POMC-Cre::Arch animal
%
% 20150602 Yuta Senzai v4: add m_minValue >
%
% INPUTS
% ch_hilus, ch_molecular: 1-based
%
% defaults
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder)
addParameter(p,'res_per',[0,inf],@isnumeric) % interval for analysis
addParameter(p,'DS2_lowThreshold',3000,@isnumeric) % 3000*0.38uV = 1.14mV
addParameter(p,'DS1_lowThreshold',3000,@isnumeric) % 3000*0.38uV = 1.14mV
addParameter(p,'DS2fil_highThreshold',3000,@isnumeric)
addParameter(p,'DS1fil_highThreshold',1500,@isnumeric)
addParameter(p,'DS1wb_highThreshold',3000,@isnumeric)
addParameter(p,'DS2_mol_threshold',500,@isnumeric)

parse(p,varargin{:})
basepath = p.Results.basepath;
res_per = p.Results.res_per; 

DS2_lowThreshold = p.Results.DS2_lowThreshold; 
DS1_lowThreshold = p.Results.DS1_lowThreshold; 
DS2fil_highThreshold = p.Results.DS2fil_highThreshold; 
DS1fil_highThreshold = p.Results.DS1fil_highThreshold; 
DS1wb_highThreshold = p.Results.DS1wb_highThreshold; 
DS2_mol_threshold = p.Results.DS2_mol_threshold; 

%% get lfp and filter
res_lfp_h = getLFP(ch_hilus,'basepath',basepath,'intervals',res_per,'noPrompts',true);
res_lfp_m = getLFP(ch_molecular,'basepath',basepath,'intervals',res_per,'noPrompts',true);

if res_lfp_m.samplingRate ~= 1250
   error('sorry, this code assumes 1250 fs') 
end

fil_res_lfp_h = bz_Filter(res_lfp_h,'passband',[2 50]);
fil_res_lfp_m = bz_Filter(res_lfp_m,'passband',[2 50]);

res_lfp_h = [res_lfp_h.timestamps,double(res_lfp_h.data)];
res_lfp_m = [res_lfp_m.timestamps,double(res_lfp_m.data)];

fil_res_lfp_h = [fil_res_lfp_h.timestamps,double(fil_res_lfp_h.data)];
fil_res_lfp_m = [fil_res_lfp_m.timestamps,double(fil_res_lfp_m.data)];

time = fil_res_lfp_h(:,1);
hm_dif = fil_res_lfp_h(:,2) - fil_res_lfp_m(:,2);


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
if DS2_start(1) > DS2_stop(1)
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
    disp(['After detection by thresholding: ',...
        num2str(length(DS2_firstPass)),' putative DS2 events.']);
end

DS2_secondPass = DS2_firstPass;


% Discard with a peak power < highThreshold for
DS2 = [];
for i = 1:size(DS2_secondPass,1)
    if DS2_secondPass(i,1)>50
        m_meanValue = mean(res_lfp_m(DS2_secondPass(i,1):DS2_secondPass(i,2),2));
        % 25/1250sec = 20ms, 40/1250 = 32 ms
        ds2criteria = mean(res_lfp_m((DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20),2))-DS2_mol_threshold; 
        
        idx = [(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20),(DS2_secondPass(i,2)+20):(DS2_secondPass(i,2)+45)];
        m_minValue =  min(res_lfp_m(idx,2));
        h_maxValue =  max(res_lfp_h(idx,2));
                                
        if (m_meanValue < ds2criteria && m_minValue > -6000 && h_maxValue < 6000)
            DS2lfp = res_lfp_h( DS2_secondPass(i,1): DS2_secondPass(i,2),:);
            [~,idx]=max(DS2lfp(:,2));
            DS2 = [DS2 ; time(DS2_secondPass(i,1)), DS2lfp(idx,1), time(DS2_secondPass(i,2))];
        end
    end
end

if isempty(DS2)
    disp('DS2 Peak thresholding failed.');
    return
else
    disp(['After peak thresholding: DS2 ', num2str(length(DS2)), ' events.']);
end


%% Detect DS1
DS1_thresholded = fil_res_lfp_h(:,2) > DS1_lowThreshold;
DS1_start = find(diff(DS1_thresholded)>0);
DS1_stop = find(diff(DS1_thresholded)<0);
% Exclude last DS if it is incomplete
if length(DS1_stop) == length(DS1_start)-1
    DS1_start = DS1_start(1:end-1);
end
% Exclude first DS if it is incomplete
if length(DS1_stop)-1 == length(DS1_start)
    DS1_stop = DS1_stop(2:end);
end
% Correct special case when both first and last DS are incomplete
if DS1_start(1) > DS1_stop(1)
    DS1_stop(1) = [];
    DS1_start(end) = [];
end
DS1_firstPass = [DS1_start,DS1_stop];
if isempty(DS1_firstPass)
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ',...
        num2str(length(DS1_firstPass)), ' putative DS1 events.']);
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
        % 25/1250sec = 20ms, 55/1250sec = 44ms
        ds1criteria = mean(res_lfp_m([(DS1_secondPass(i,1)-80):(DS1_secondPass(i,1)-55) ],2))+500; 
        m_maxValue = max(res_lfp_m([DS1_secondPass(i,1)-55:DS1_secondPass(i,2)+55],2));
        m_minValue = min(res_lfp_m([DS1_secondPass(i,1)-55:DS1_secondPass(i,2)+55],2));
        % make sure period does not exceed the size of res_lfp
        if size(res_lfp_h,1)>=(DS1_secondPass(i,2)+125)
            % 100/1250sec = 80ms, 125/1250sec = 100ms
            h_postValue = mean(res_lfp_h([(DS1_secondPass(i,1)+100):(DS1_secondPass(i,2)+125)],2)); 
            % 25/1250sec = 20ms, 55/1250sec = 44ms
            ds1postcriteria = mean(res_lfp_h([(DS1_secondPass(i,1)-80):(DS1_secondPass(i,1)-55) ],2))+2000; 
            if (df_maxValue > DS1fil_highThreshold) && (h_maxValue > DS1wb_highThreshold) && (m_meanValue > ds1criteria) ...
                    && (h_postValue < ds1postcriteria) && (m_minValue > -3000) && (m_maxValue < 3000)
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

% .event.mat (buzcode & cell_explorer standard)
clear DS1 DS2

DS1.timestamps = DS1triad(:,[1 3]);
DS1.peaks = DS1triad(:,2);
DS1peakIDCs = round(DS1triad(:,2)*1250);
% .00038 = conversion to mV
DS1.amplitudes = (res_lfp_h(DS1peakIDCs,2)-res_lfp_m(DS1peakIDCs,2))*.00038; 
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
DS1.duration = DS1triad(:,3) - DS1triad(:,1);
DS1.center = DS1triad(:,1) + DS1triad(:,2) - DS1triad(:,1);

DS1.detectorinfo.detectorname = 'DetectDSpikes_v4';
DS1.detectorinfo.detectionparms = [];
DS1.detectorinfo.detectionintervals = res_per;
DS1.detectorinfo.detectiondate = datetime('today');
DS1.detectorinfo.ml_channel = ch_molecular;
DS1.detectorinfo.h_channel = ch_hilus;

DS2.timestamps = DS2triad(:,[1 3]);
DS2.peaks = DS2triad(:,2);
DS2peakIDCs = round(DS2triad(:,2)*1250);
% .00038 = conversion to mV
DS2.amplitudes = (res_lfp_h(DS2peakIDCs,2)-res_lfp_m(DS2peakIDCs,2))*.00038; 
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
DS2.duration = DS2triad(:,3) - DS2triad(:,1);
DS2.center = DS2triad(:,1) + DS2triad(:,2) - DS2triad(:,1);

DS2.detectorinfo.detectorname = 'DetectDSpikes_v4';
DS2.detectorinfo.detectionparms = [];
DS2.detectorinfo.detectionintervals = res_per;
DS2.detectorinfo.detectiondate = datetime('today');
DS2.detectorinfo.ml_channel = ch_molecular;
DS2.detectorinfo.h_channel = ch_hilus;

basename = basenameFromBasepath(basepath);
save(fullfile(basepath,[basename '.DS1.events.mat']),'DS1');
save(fullfile(basepath,[basename '.DS2.events.mat']),'DS2');
end