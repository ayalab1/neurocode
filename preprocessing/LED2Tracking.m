
function [tracking] = LED2Tracking(aviFile,varargin)
% Get LED tracking
%
% USAGE
%
%   [behavior] = LED2Tracking(varagin)
%
% INPUTS
%   aviFile     Avi format video. If not provided, look for it in the
%   basepath f
%   It requires an avi format video and a digitalin.dat file in the
%   basepath folder.
% 
%   (OPTIONAL)
%   basePath       -(default: pwd) basePath for the recording file, in
%                   buzcode format.
%   roiTracking    - 2 x R, where 1C is x and 2C is y. By default it
%                   considers the whole video. With the option 'manual' allows to draw
%                   a ROI.
%   fs             - video sampling rate, default = 30 Hz
%   roiLED         - 2 x R, where 1C is x and 2C is y.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveFrames     - Creates mat file containin all frames (default false).
%   forceReload    - default false.
%   RGBChannel     - 'r', 'g', 'b', or any combination (ej. 'rgb'); default
%                    'r' TO DO!!
%   bazlerTTL      - Rx1 with timestamps from bazler ttl pulses. If not
%                   provided try to extract ttl pulses from digitalIn
%                   channel 1. If it fails, gives video time.
%   saveMat        - default true
%   artifactThreshold - max allow movements per frame (in cm, default 3).
%                   Disabled if not convFact is provided.
% 
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   x               - x position in cm/ normalize
%   y               - y position in cm/ normalize
%   timestamps      - in seconds, if Bazler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%   samplingRate    - in Hz
%   averageFrame    - 
%   
%
%   aza oliva, 2021, 
%   based on previous one from Manu Valero 2019, and  Adrien Peyrache 2015
%   adapted to neurocode and setups in ayalab (Nov - 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'fs',30,@isnumeric);
addParameter(p,'artifactThreshold',10,@isnumeric);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiTracking',[],@ismatrix);
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveFrames',true,@islogical)
addParameter(p,'verbose',false,@islogical);
addParameter(p,'thresh',.75,@isnumeric) %0.98
addParameter(p,'bazlerTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'pixelSizeEstimate',[1 30],@ismatrix);
% addParameter(p,'RGBChannel',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
artifactThreshold = p.Results.artifactThreshold;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
forceReload = p.Results.forceReload;
saveFrames = p.Results.saveFrames;
verbose = p.Results.verbose;
thresh = p.Results.thresh;
bazlerTtl = p.Results.bazlerTTL;
saveMat = p.Results.saveMat;
pixelSizeEstimate = p.Results.pixelSizeEstimate;

% RGBChannel = p.Results.RGBChannel;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) || forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

if ~exist('aviFile') || isempty(aviFile)
    if ~isempty(dir([basepath filesep '*Basler*avi']))
        aviFile = dir([basepath filesep '*Basler*avi']); 
        aviFile = erase(aviFile.name,'.avi');
    else
        warning('No video file!!');
        tracking = [];
        return
    end
end
% attention, for now it only loads the red channel from the video!!
if ~exist([basepath filesep aviFile '.mat'],'file') && ~forceReload
    disp('Get average frame...');
    videoObj = VideoReader([aviFile '.avi']);   % get video
    numFrames = get(videoObj, 'NumFrames');
    clear temp
    batches = 1:2000:numFrames;
    batches = [batches numFrames+1];
    frames.r = [];frames.b = [];
    tic
    f = waitbar(0,'Getting frames...');
    for ii = 1:length(batches)-1
        waitbar(ii/length(batches),f)
        temp_frames = read(videoObj,[batches(ii) batches(ii+1)-1]);        % get all frames
        frames.r = cat(3,frames.r,squeeze(temp_frames(:,:,1,:)));          % get red shiny pixels
        frames.b = cat(3,frames.b,squeeze(temp_frames(:,:,2,:)));          % get red shiny pixels
    end
    close(f)
    toc
    
    if saveFrames
        disp('Saving frames...');
        save([basepath filesep aviFile '.mat'],'frames','-v7.3');
    end
else
    disp('Loading frames from mat file...');
    load([basepath filesep aviFile '.mat'],'frames');
end

% get average frame - plot is shity some times, usually not a problem
average_frame = mean(frames.r,3);                                          

% deal with the ROI for the tracking 
cd(basepath); cd ..; upBasepath = pwd; cd(basepath);
if exist([basepath filesep 'roiTracking.mat'],'file')
    load([basepath filesep 'roiTracking.mat'],'roiTracking');
elseif exist([upBasepath filesep 'roiTracking.mat'],'file')
    load([upBasepath filesep 'roiTracking.mat'],'roiTracking');
    disp('ROI tracking from master folder... copying locally...');
    save([basepath filesep 'roiTracking.mat'],'roiTracking');
elseif isempty(roiTracking)
    roiTracking = [1 1 size(frames.r,2) size(frames.r,2) 1 ;1 size(frames.r,1) size(frames.r,1) 1 1 ]';
elseif ischar(roiTracking) && strcmpi(roiTracking,'manual')
    disp('Draw ROI for tracking...');
    h1 = figure;
    imshow(average_frame);
    title('Draw ROI for tracking... (include all posible locations)','FontWeight','normal');
    roi = drawpolygon;
    roiTracking = [roi.Position; roi.Position(1,:)];
    save([basepath filesep 'roiTracking.mat'],'roiTracking');
    close(h1);
end

% deal with the ROI for the LED
if exist([basepath filesep 'roiLED.mat'],'file')
    load([basepath filesep 'roiLED.mat'],'roiLED');
elseif exist([upBasepath filesep 'roiLED.mat'],'file')
    load([upBasepath filesep 'roiLED.mat'],'roiLED');
    disp('ROI LED from master folder... copying locally...');
    save([basepath filesep 'roiLED.mat'],'roiLED');
elseif ischar(roiLED) && strcmpi(roiLED,'manual')
    disp('Draw ROI for LED...');
    h1 = figure;
    imshow(average_frame);
    roi = drawpolygon;
    roiLED = [roi.Position; roi.Position(1,:)];
    save([basepath filesep 'roiLED.mat'],'roiLED');
    close(h1);
end

% if convFact not provided, normalize to 1 along the longest axis
if isempty(convFact)                            
    convFact = 1/max([size(frames.r,1) size(frames.r,2)]);
    artifactThreshold = Inf;
end
xMaze = [0 size(frames.r,2) * convFact];
yMaze = [0 size(frames.r,1) * convFact];

% save ROI figure
h1 = figure;
hold on
imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 4*mean(average_frame(:))]);
set(gca,'YDir','normal', 'TickDir','out');
p = plot(roiTracking(:,1)*convFact, roiTracking(:,2)*convFact,'r','LineWidth',2);
axis tight
legend(p,'Tracking ROI');
xlabel('Normalize/ cm');
mkdir('Behavior');
saveas(h1,'Behavior\MazeROI.png');
if ~verbose
    close(h1);
end

%% DETECT LED POSITION
bw = uint8(poly2mask(roiTracking(:,1),roiTracking(:,2),size(frames.r,1),size(frames.r,2)));
disp('Detect LED position...');
thr_fr = thresh * 255;
if ~verbose
    tic
    f = waitbar(0,'Detecting LED position...');
    for ii = 1:size(frames.r,3)
        waitbar(ii/size(frames.r,3),f)
        
        %red LED
        fr = frames.r(:,:,ii).*bw;
        bin_fr = imbinarize(double(fr),thr_fr); %
        bin_fr = bwareafilt(bin_fr,pixelSizeEstimate);
        stats_fr = regionprops(bin_fr);
        maxBlobr = find([stats_fr.Area]== max([stats_fr.Area]),1);
        
        %blue LED
        fb = frames.b(:,:,ii).*bw;
        bin_fb = imbinarize(double(fb),thr_fr); %
        bin_fb = bwareafilt(bin_fb,pixelSizeEstimate);
        stats_fb = regionprops(bin_fb);
        maxBlobb = find([stats_fb.Area]== max([stats_fb.Area]),1);
        
        %red LED
        if ~isempty(maxBlobr)
            sz_fr(ii) = stats_fr(maxBlobr).Area;
            Rr_x(ii) = stats_fr(maxBlobr).Centroid(1);
            Rr_y(ii) = stats_fr(maxBlobr).Centroid(2);
        else
            sz_fr(ii) = NaN;
            Rr_x(ii) = NaN;
            Rr_y(ii) = NaN;
        end
        
        %blue LED
        if ~isempty(maxBlobb)
            sz_fb(ii) = stats_fb(maxBlobb).Area;
            Rb_x(ii) = stats_fb(maxBlobb).Centroid(1);
            Rb_y(ii) = stats_fb(maxBlobb).Centroid(2);
        else
            sz_fb(ii) = NaN;
            Rb_x(ii) = NaN;
            Rb_y(ii) = NaN;
        end
        
    end
    close(f)
    toc
else
    h1 = figure;
    hold on
    tic
    for ii = 1:size(frames.r,3)
        
        %red LED
        fr = frames.r(:,:,ii).*bw;
        bin_fr = imbinarize(double(fr),thr_fr); %
        bin_fr = bwareafilt(bin_fr,pixelSizeEstimate);
        stats_fr = regionprops(bin_fr);
        maxBlobr = find([stats_fr.Area]== max([stats_fr.Area]),1);
        
        %blue LED
        fb = frames.r(:,:,ii).*bw;
        bin_fb = imbinarize(double(fb),thr_fr); %
        bin_fb = bwareafilt(bin_fb,pixelSizeEstimate);
        stats_fb = regionprops(bin_fb);
        maxBlobb = find([stats_fb.Area]== max([stats_fb.Area]),1);
        
        %red LED
        if ~isempty(maxBlobr)
            sz_fr(ii) = stats_fr(maxBlobr).Area;
            Rr_x(ii) = stats_fr(maxBlobr).Centroid(1);
            Rr_y(ii) = stats_fr(maxBlobr).Centroid(2);
            cla
            imagesc(fr)
            plot(Rr_x(ii),Rr_y(ii),'or')
            drawnow;
        else
            sz_fr(ii) = NaN;
            Rr_x(ii) = NaN;
            Rr_y(ii) = NaN;
        end
        
        %blue LED
        if ~isempty(maxBlobb)
            sz_fb(ii) = stats_fb(maxBlobb).Area;
            Rb_x(ii) = stats_fb(maxBlobb).Centroid(1);
            Rb_y(ii) = stats_fb(maxBlobb).Centroid(2);
            cla
            imagesc(fb)
            plot(Rb_x(ii),Rb_y(ii),'or')
            drawnow;
        else
            sz_fb(ii) = NaN;
            Rb_x(ii) = NaN;
            Rb_y(ii) = NaN;
        end
        
    end
    toc
    close(h1);
end

pos_r = [Rr_x; Rr_y]';
pos_b = [Rb_x; Rb_y]';

if isempty(find(~isnan(pos_r), 1)) && isempty(find(~isnan(pos_b), 1))
    warning('Not finding positions... usually a problem of proper threshold to detect LED');
end

%% postprocessing of LED position

pos1 = pos_r * convFact;                                  % cm or normalized
pos2 = pos_b * convFact;                                  % cm or normalized
art = find(sum(abs(diff(pos1))>artifactThreshold,2))+1;  % remove artefacs as movement > 10cm/frame
pos1(art,:) = NaN;
art2 = find(sum(abs(diff(pos2))>artifactThreshold,2))+1;  % remove artefacs as movement > 10cm/frame
pos2(art,:) = NaN;

xt = linspace(0,size(pos1,1)/fs,size(pos1,1));            
% kalman filter is not working - NEED TO BE FIXED
% [t,x,y,vx,vy,ax,ay] = trajectory_kalman_filter(pos(:,1)',pos(:,2)',xt,0);
x1=pos1(:,1);y1=pos1(:,2); t=xt;
x2=pos2(:,1);y2=pos1(:,2);
art = find(sum(abs(diff([x1 y1]))>artifactThreshold,2))+1;
art = find(sum(abs(diff([x2 y2]))>artifactThreshold,2))+1;
art = [art - 2 art - 1 art art + 1 art + 2];
art2 = [art2 - 2 art2 - 1 art2 art2 + 1 art2 + 2];
x1(art(:)) = NaN; y1(art(:)) = NaN;
x2(art2(:)) = NaN; y2(art2(:)) = NaN;
F = fillmissing([x1 y1],'linear');
F2 = fillmissing([x2 y2],'linear');
x1 = F(:,1); y1 = F(:,2);
x2 = F2(:,1); y2 = F2(:,2);

h2 = figure;
hold on
imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 .7]);
freezeColors;
scatter(x1,y1,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
caxis([t(1) t(end)])
xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
xlim(xMaze); ylim(yMaze);
mkdir('Behavior');
saveas(h2,'Behavior\trajectory.png');
if ~verbose
    close(h2);
end

%% detect LED pulses for sync
if ~isempty(roiLED)
    disp('Detect LED for sync...');
    bwLED = uint8(poly2mask(roiLED(:,1),roiLED(:,2),size(frames,1),size(frames,2)));
    parfor ii = 1:size(frames,4)
        fr = double(frames(:,:,1,ii).*bwLED);
        fr(fr==0) = NaN;
        sync(ii) = nansum(fr(:)); 
    end

    sync = sync.^2;
    syncBin = (sync>mean(sync)); % binarize signal
    locsA = find(diff(syncBin)==1)/fs; % start of pulses
    locsB = find(diff(syncBin)==-1)/fs; % end of pulses
    pul = locsA(1:min([length(locsA) length(locsB)]));
    for ii = 1 : size(pul,2) % pair begining and end of the pulse
        if sum(locsB > pul(1,ii)) > 0
            pul(2,ii) =  locsB(find(locsB - pul(1,ii) ==...
                min(locsB(locsB > pul(1,ii)) - pul(1,ii))));
        else
            pul(2,ii) = nan;
        end
    end
else
    sync = []; pul = [];
end

%% Get basler TTL
% digital_in legend: 1. Basler, 2. maze LEd, 3. Left Alternation, 4.Righ
% Alternation, 5. Home Delay, 6. Is alternation forzed?
if isempty(bazlerTtl)
    digitalIn = getDigitalIn;
    bazlerTtl = digitalIn.timestampsOn{1};
end

% match basler frames con ttl pulses
if length(bazlerTtl) == length(x1)
    disp('Number of frames match!!');
elseif length(bazlerTtl) > length(x1) && length(bazlerTtl) <= length(x1) + 15 * 1 
    fprintf('%3.i frames were dropped, probably at the end of the recording. Skipping... \n',...
        length(bazlerTtl) - length(x1));
    bazlerTtl = bazlerTtl(1:length(x1));
elseif length(bazlerTtl) < length(x1) && (length(x1)-length(bazlerTtl)) < 30 * 4
    fprintf('%3.i video frames without TTL... Was the recording switched off before the camera?. Skipping... \n',...
        length(x1) - length(bazlerTtl));
    x1 = x1(1:length(bazlerTtl));
    y1 = y1(1:length(bazlerTtl));
    x2 = x2(1:length(bazlerTtl));
    y2 = y2(1:length(bazlerTtl));
    vx = vx(1:length(bazlerTtl));
    vy = vy(1:length(bazlerTtl));
    ax = ax(1:length(bazlerTtl));
    ay = ay(1:length(bazlerTtl)); 
elseif isempty(bazlerTtl)
    bazlerTtl = xt;
elseif abs(length(x1)-length(bazlerTtl)) > 15 * 1 && size(digitalIn.timestampsOn,2)> 4
    fprintf('%3.i frames were dropped, possibly at the beginning of the recording. Aligning timestamps to the first IR TTL... \n',...
        length(bazlerTtl) - length(x1));
    f1 = figure;
    hold on
    imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 4*mean(average_frame(:))]); axis tight
    set(gca,'Ydir','reverse');    
    bazlerTtl = bazlerTtl(1:length(x));
    
else
%    keyboard;
    error('Frames do not match for more than 5 seconds!! Trying to sync LED pulses');
    if length(digitalIn.timestampsOn{2}) == length(sync_signal)
         disp('Using sync LED pulses...');
         keyboard; % to do!!!
    end
end

[~,fbasename,~]=fileparts(pwd);

tracking.position.x1 =x1;
tracking.position.y1 = y1;
tracking.position.x2 =x2;
tracking.position.y2 = y2;
tracking.position.z = [];
tracking.description = '';
tracking.timestamps = bazlerTtl;
tracking.originalTimestamps = [];
tracking.folder = fbasename;
tracking.sync.sync = sync;
tracking.sync.timestamps = pul;
tracking.samplingRate = fs;
tracking.avFrame.r = average_frame;
tracking.avFrame.xSize = xMaze;
tracking.avFrame.ySize = yMaze;
tracking.roi.roiTracking = roiTracking;
tracking.roi.roiLED = roiLED;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

% maze is 48 x 67 (exterior wall to exterior wall). Virtual maze is 580 x
% 420 pixels. Conv factor is ~ 0.1143 - 0.1155 cm/pix 
