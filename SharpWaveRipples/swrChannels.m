
function rippleChannels = swrChannels(varargin)
%swrChannels - optimize ripple, SW, and noise channel selection
% rippleChannel = swrChannels(varargin)

% Find ripple channel of interest (ripple, SWP and noise channels)

% INPUT
%   <options>       optional list of property-value pairs (see table below)
%   basepath        Basepath containing...
%   discardShanks   Default [].
%   probesNumber    Default 1.
%   noPrompts       Default, true
%   saveMat         Default, true.
%   force           Default, false
%
%   F.Sharif 2020. Converted to buzcode from Channel_Info by MV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'discardShanks',[],@isnumeric);
addParameter(p,'probesNumber',1,@isnumeric);
addParameter(p,'noPrompts',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
parse(p,varargin{:});
basepath = p.Results.basepath;
discardShanks = p.Results.discardShanks;
probesNumber = p.Results.probesNumber;
noPrompts = p.Results.noPrompts;
saveMat = p.Results.saveMat;
force = p.Results.force;

prevPath = pwd;
cd(basepath);

filename = dir('*.channelinfo.ripples.mat');
if ~isempty(filename) && ~force
    disp('Ripple channel already detected! Loading file.');
    load(filename.name);
    return
end

% main code (Channel_Info)
basename = basenameFromBasepath(basepath);
load([basename '.session.mat']);
Anatomical_groups=session.extracellular.electrodeGroups.channels;
disp('Reading .XML /  Anatomical Groups ')
disp('...                                  ')


SHANKS=[];
for GroupsNumber=1:size(Anatomical_groups,2)
    SHANKS{GroupsNumber}= cell2mat(Anatomical_groups(1, GroupsNumber));
    disp(['Group ' num2str(GroupsNumber) ' = '  '[' num2str(SHANKS{1,GroupsNumber}) ']' ]);
end

channelnumber_correction = 0;

if length(SHANKS)> 9 | probesNumber > 1
    disp('...                                  ')
    disp('There are more than 1 prob in this recording session')
    Seperate_Probs=input('Do you like to find individual refrence channel for each prob?  (if yes press 1 else press 0) ');
    if Seperate_Probs==1
        GroupsNumber_User=input('insert .XML groups numbers (Shank numbers) belong to prob1 .ex [1:8] or [ 1 2 3 ...]')
        CorticalChannel_included=input('Are cortical cahnnels included in the above groups (if yes press 1 else press 0) ');
        if CorticalChannel_included==0
            CorticalChannel_User=input('What are cortical channels needed to be merged to this prob .ex [128 129 ...] ')
        else
            CorticalChannel_User=[];
        end

    end
else
    Seperate_Probs=-1;
    disp('                                ')
    disp('1 probs is detectted')
    disp('                                ')
    disp('                                ')
end
    
SHANKS_PerProb=[];
if Seperate_Probs==1
    for i=1:length(GroupsNumber_User)
        SHANKS_PerProb{1,i}= SHANKS{1,GroupsNumber_User(i)}+channelnumber_correction;
        disp(['NewGroup ' num2str(i) ' = '  '[' num2str(SHANKS_PerProb{1,i}) ']' ])
    end
    if isempty(CorticalChannel_User)==0
        SHANKS_PerProb{1,i+1}=CorticalChannel_User+channelnumber_correction;
        disp(['NewGroup ' num2str(i+1) ' = '  '[' num2str(SHANKS_PerProb{1,i+1}) ']' ])
    end
else
    for i=1:length(SHANKS)
        SHANKS_PerProb{1,i}= SHANKS{1,i}+channelnumber_correction;
        disp(['NewGroup ' num2str(i) ' = '  '[' num2str(SHANKS_PerProb{1,i}) ']' ])
    end
end
clear SHANKS
SHANKS=SHANKS_PerProb;

for ii = 1:length(discardShanks)
    SHANKS(discardShanks) = [];
end

%% Calculate Signal to noise ratio in all of the channels ##################

Ripple_window_start=120;
Ripple_window_end=180;
[lfp] = getLFP('all');
[bbb aaa]=butter(4,[Ripple_window_start Ripple_window_end]/lfp.samplingRate*2,'bandpass');

Refrence_chan=[];
% Not Necessary, only for test ############################################
Signal2Noise=[];
StandDev=[];
% #########################################################################

    for shk=1:length(SHANKS)

        clear var ChannelsPerShank SNR

        for CH=1:length(SHANKS{shk})

            filt = filtfilt(bbb,aaa,double(lfp.data(:,SHANKS{shk}(CH))));
            Amp = fastrms(filt,15);
            % Calculate signal to noise ratio #################################
            Amp_Sort=sort(Amp);
            percen=(1e-3)/100;
            Alpha=floor(length(Amp_Sort)*percen);
            minAmp=mean(Amp_Sort(1:Alpha));
            maxAmp=mean(Amp_Sort(end-Alpha:end));
            SNR(CH) = 10.*log((maxAmp-minAmp)./std(Amp));

            % Not Necessary, only for test ####################################
            Signal2Noise=[Signal2Noise;SNR(CH)];
            StandDev=[StandDev;1./std(Amp)];
            % #################################################################
        end

        SNR=SNR./max(SNR);
        [minVals, locs] = min(SNR);
        Refrence_chan(shk,:) = [minVals locs SHANKS{shk}(locs)];

    end

[~,RefrenceShnak]=min(Refrence_chan(:,1));
RefrenceRippleChannel_test=Refrence_chan(RefrenceShnak,3);
disp(['Reference Ripple Channel for test is : ' num2str(RefrenceRippleChannel_test) ]);

%% Get Ripple
[ripples] = FindRipples('basepath',basepath,'channel',RefrenceRippleChannel_test);
Win=70;
LfpSamplingrate = lfp.samplingRate;
% Removing short startting and the end ripples
ripples.peaks = ripples.peaks(ripples.peaks*LfpSamplingrate>Win+1 & ripples.peaks*LfpSamplingrate<length(lfp.timestamps)-Win+1);

%% Calculate Ripple power ##################################################

Rippl_Matrix=[];
Ripples_Power=[];
Ripples_Power_Matrix=[];
Ripples_CSD=[];

for shk=1:length(SHANKS)

    for CH=1:length(SHANKS{shk})

        clear var All_Ripple_Avg  ripple_ave Power

        eeg=single(lfp.data(:,SHANKS{1, shk}(CH)));

        for i = 1:size(ripples.peaks,1)
            ripple_ave(i,:) = eeg(round(ripples.peaks(i)*LfpSamplingrate)-Win:round(ripples.peaks(i)*LfpSamplingrate)+Win,:);
        end

        All_Ripple_Avg=double(mean(ripple_ave));

        Frq=120:180;
        scale=frq2scal(Frq,LfpSamplingrate);
        S=cwt(All_Ripple_Avg,scale,'morl');
        g_baseline= (envelop(S.*S))';
        Power=mean(g_baseline(50:100,:),2);

        Rippl_Matrix{1,shk}(CH,:)=All_Ripple_Avg;
        Ripples_Power=[Ripples_Power; mean(Power) SHANKS{1,shk}(CH)] ;
        Ripples_Power_Matrix{1,shk}(CH,:)=Power;
        Ripples_CSD{1,shk}(CH,:)=[sum(diff(smooth1D(All_Ripple_Avg(50:70),10),2)) sum(diff(smooth1D(All_Ripple_Avg(70:100),10),2)) SHANKS{1,shk}(CH) ];
    end

end

%% find Deep and superficial channels ######################################

Deep_Sup=[];
con_sum_all=[];
con_direction_all=[];

for shk=1:size(SHANKS,2)
    clear Reversal_channel
    con_direction=Ripples_CSD{1,shk}(:,1).*Ripples_CSD{1,shk}(:,2);
    con_sum=Ripples_CSD{1,shk}(:,1)+Ripples_CSD{1,shk}(:,2);
    con_sum_all=[con_sum_all;con_sum Ripples_CSD{1,shk}(:,3)];
    con_direction_all=[con_direction_all; con_direction Ripples_CSD{1,shk}(:,3)];
    nd=find(con_direction<0);
    ndx=find(con_sum<0);

    if isempty(nd)==0
        Reversal_channel=nd(end);
        Deep_Sup{1,shk}(1:Reversal_channel,1)=Reversal_channel-(1:Reversal_channel);
        if Reversal_channel< size(SHANKS{1,shk},2)
            Deep_Sup{1,shk}(Reversal_channel+1:size(SHANKS{1,shk},2),1)=Reversal_channel-(Reversal_channel+1:size(SHANKS{1,shk},2));
        end       
    elseif isempty(ndx)==0
        Reversal_channel=ndx(end);
        Deep_Sup{1,shk}(1:Reversal_channel,1)=Reversal_channel-(1:Reversal_channel);
        if Reversal_channel< size(SHANKS{1,shk},2)
            Deep_Sup{1,shk}(Reversal_channel+1:size(SHANKS{1,shk},2),1)=Reversal_channel-(Reversal_channel+1:size(SHANKS{1,shk},2));
        end
    else
        Deep_Sup{1,shk}(:,1)=-(1:size(SHANKS{1,shk},2));
    end

end

%%
% Assigning  ripple , SWR , and noise channel #############################
[~,nd_Ripple]=max(Ripples_Power(:,1));
Rip_chnl=Ripples_Power(nd_Ripple,2);

con_direction_all(con_direction_all(:,1)<0,1)=NaN;
[~,nd_noise]=min(con_direction_all(:,1));
noise_chnl=con_direction_all(nd_noise,2);

[~,nd_SWR]=max(con_sum_all(:,1));
SWR_chnl=con_sum_all(nd_SWR,2);

%% 1-Plot ripple layout

figure('position',[200 115 1300 800])
subplot(4,1,[1:3])
Sh_spacing=200;
chan_spacing=800;
ylabel_p=[];
K=0;
for shk=1:length(SHANKS)
    
    clear var SD
    
    for CH=1:length(SHANKS{shk})
        K=K+1;
        
        clear var Ripple_channel Cr
        
        Ripple_Y=Rippl_Matrix{shk}(CH,:)-(CH-1)*chan_spacing;
        Ripple_X=[1:length(Ripple_Y)]+(shk-1)*Sh_spacing;
        
        SD=Deep_Sup{shk}(CH);
     
        % Deep sup layer Colors ###########################################
        if SD <0
            Cr=[0 146 146]./256;
        else
            Cr=[0 109 219]./256;
        end
        plot(Ripple_X,Ripple_Y,'color',Cr,'linewidth',1);
        hold on
        % Channel Colors ##################################################
        if SHANKS{1,shk}(CH)==Rip_chnl
            channelname='Ripple';
            plot(Ripple_X,Ripple_Y,'color','r','linewidth',1,'LineStyle','-.');
            text(Ripple_X(end)-20,Ripple_Y(1)+2*Sh_spacing,channelname,'color','r','fontsize',10)
        elseif SHANKS{1,shk}(CH)==SWR_chnl
            channelname='SWR';
            plot(Ripple_X,Ripple_Y,'color','r','linewidth',1,'LineStyle','-.');
            text(Ripple_X(end),Ripple_Y(1)+2*Sh_spacing,channelname,'color','r','fontsize',10)
        elseif SHANKS{1,shk}(CH)==noise_chnl
            channelname='noise';
            plot(Ripple_X,Ripple_Y,'color','r','linewidth',1,'LineStyle','-.');
            text(Ripple_X(end),Ripple_Y(1)+2*Sh_spacing,channelname,'color','r','fontsize',10)
        end

        %Type channel number and deep sup
        %##################################################################
        if SD <0
        text(Ripple_X(1)-30,Ripple_Y(1)+2*Sh_spacing,['Sup' num2str([SHANKS{1,shk}(CH)])])
        else
        text(Ripple_X(1)-30,Ripple_Y(1)+2*Sh_spacing,['Deep' num2str([SHANKS{1,shk}(CH)])])
        end
        %Type test number #################################################
        hold on
    end
end

set(gca,'YTick',[])
set(gca,'visible','off')
set(gca,'color','w')

%% 2-Plot Signal to noise ratio in all of the channels
% figure('position',[200 700 500 270])
subplot(4,3,10)
Sh_spacing=1;
chan_spacing=.5;
yticks_p=[];

for shk=1:length(SHANKS)
    
    clear var ChannelsPerShank
    
    for CH=1:length(SHANKS{shk})
        
        if SHANKS{shk}(CH)==RefrenceRippleChannel_test
            Cr='m';
        elseif CH==Refrence_chan(shk,2)
            Cr='r';
        else
            Cr=[222, 222, 222]/256;
        end
        
        %Plot channels ####################################################
        recX=shk*2;
        recY=CH-chan_spacing/2;
        Ripple_X=[recX recX+Sh_spacing recX+Sh_spacing recX];
        Ripple_Y=[recY recY recY+chan_spacing  recY+chan_spacing];
        fill(Ripple_X,Ripple_Y,Cr,'edgecolor','none');
        hold on
    end
    
    % x and y ticks label position ########################################
    xticks_p(shk)=recX+Sh_spacing/2;
    yticks_p=[yticks_p;length(SHANKS{shk})];
    % import channel number having minimum NSR at each shank ##############
    MinSNRchannel=SHANKS{1,shk}(Refrence_chan(shk,2));
    text(recX,0,num2str([MinSNRchannel]))
    
end
xticks(xticks_p)
xticklabels(round(Refrence_chan(:,1),4))
xlabel('Minimum SNR per shank','fontsize', 14)
yticks(1:max(yticks_p))
ylim([0 max(yticks_p)+1])
ylabel('Channel Number','fontsize', 14)
axis ij
title('XML-Anatomical Groups Number','fontsize', 10)
box('off')

%% 3-Plot Ripple-noise correlation

a=[];
a=Ripples_Power(:,1);
a=a./max(a);
b=[];
b=Signal2Noise(:,1);
b=b./max(b);
[r,p]=corrcoef(a,b);
R=r(1,2);
P=p(1,2);
    
% figure('position',[700 700 300 270])
subplot(4,3,11)
plot(a,b,'ok','linewidth',3)
hold on;
h = lsline;
set(h,'color','r','linewidth',2)
xlabel('Ripple power')
ylabel('SNR')
title (['r = ' num2str(R) ['    pval = ' num2str(P)]] )

%% 4-plot Power layout
colmn=size(Ripples_Power_Matrix,2);
for i=1:colmn
    ChannelNumber(i)=size(Ripples_Power_Matrix{1, i},1);
end
rows=max(ChannelNumber);
Ripple_length=size(Ripples_Power_Matrix{1, 1},2);
Power_Matrix=zeros(rows,colmn*Ripple_length);

for shk=1:colmn
    begin_ndx=(shk-1)*Ripple_length+1;
    end_ndx=shk*Ripple_length;
    rows_number=ChannelNumber(shk);
    Power_Matrix(1:rows_number,begin_ndx:end_ndx)= Ripples_Power_Matrix{1, shk};
end

% figure('position',[1000 700 500 270])
subplot(4,3,12)
imagesc(Power_Matrix)
text(1,1,'1','color','w')
colormap jet
xlabel('Shanks')
ylabel('Channels')
set(gca,'YTick',[])
set(gca,'XTick',[])

if ~noPrompts
    Bout_2 = input(sprintf('Do you like the selected channels?  (if yes press 1 else press any key)   '));
    if Bout_2==1
        disp(['Ripple_Channel=' num2str(Rip_chnl) ' Sharpwave_Channel=' num2str(SWR_chnl) '  Noise_Chnnel='  num2str(noise_chnl) ])
        disp('Done!')
    else

        Bout_3 = input(sprintf('Insert Ripple,Swr,and noise as : [Ripple_Chnnel# Swr_Chnnel# noise_Chnnel#]'));
        disp(['Ripple_Channel=' num2str(Bout_3(1)) ' Sharpwave_Channel=' num2str(Bout_3(2)) '  noise_Chnnel=' num2str(Bout_3(3)) ])
    end
end

rippleChannels.Ripple_Channel=Rip_chnl;
rippleChannels.Sharpwave_Channel=SWR_chnl;
rippleChannels.Noise_Channel=noise_chnl;
rippleChannels.Deep_Sup=Deep_Sup;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.channelinfo.ripples.mat'],'rippleChannels');
end
    

cd(prevPath);
end

%% Additional functions

function [maxsind, maxsvalues] = LocalMaxima(x)

    nPoints = length(x);
    Middle = x(2:(nPoints-1));
    Left = x(1:(nPoints-2));
    Right = x(3:nPoints);
    maxsind = 1+find(Middle > Left & Middle > Right);
    maxsvalues = x(maxsind);
    
end

function scal = frq2scal(freq,samplingrate)

    S = 2:0.001:200;
    sample_period = 1/samplingrate;
    f = scal2frq(S,'morl',sample_period);
    for ii = 1:length(freq)
        [c,ndx]=min(abs(f-freq(ii)));
        scal(ii) = S(ndx);
    end
    
end

function Senv = envelop(S)

    npts = length(S(1,:));
    for ii = 1:length(S(:,1))
        [maxI,maxV] = LocalMaxima(S(ii,:));
        if length(maxI)>1
            Senv(ii,:) = interp1([1 maxI npts],[S(ii,1) maxV S(ii,end)],1:npts,'spline');
        else
            Senv(ii,:) = S(ii,:);
        end
    end
    
end

function Sdata=smooth1D(data,halfwidth,dim)

    Fullwin=(-halfwidth:halfwidth)';
    Svector=exp(-Fullwin.^2/(halfwidth/2)^2);

    [a,b]=size(data);
    if a>1&b>1
        if dim==1
            for ii=1:length(data(1,:))
            smoothed=conv(data(:,ii),Svector)/sum(Svector);
            Sdata(:,ii)=smoothed((halfwidth+1):(end-halfwidth));
            end
        else
            for ii=1:length(data(:,1))
            smoothed=conv(data(ii,:),Svector)/sum(Svector);
            Sdata(ii,:)=smoothed((halfwidth+1):(end-halfwidth));
            end
        end
    else

        smoothed=conv(data,Svector)/sum(Svector);
        Sdata=smoothed((halfwidth+1):(end-halfwidth));

    end
    
end