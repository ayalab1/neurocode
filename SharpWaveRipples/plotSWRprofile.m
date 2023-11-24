
function  plotSWRprofile

basename = basenameFromBasepath(pwd);
load([basename '.session.mat']);
Anatomical_groups=session.extracellular.electrodeGroups.channels;

SHANKS=[];
for GroupsNumber=1:size(Anatomical_groups,2)
    SHANKS{GroupsNumber}= cell2mat(Anatomical_groups(1, GroupsNumber));
    disp(['Group ' num2str(GroupsNumber) ' = '  '[' num2str(SHANKS{1,GroupsNumber}) ']' ]);
end

load([basename '.ripples.events.mat']);
Win=70;LfpSamplingrate = 1250;
ripples.peaks = ripples.peaks(ripples.peaks*LfpSamplingrate>Win+1 & ripples.peaks*LfpSamplingrate-Win+1);

[lfp] = getLFP('all');
for shk=1:length(SHANKS)
    for CH=1:length(SHANKS{shk})
        eeg=single(lfp.data(:,SHANKS{1, shk}(CH)));
        for i = 1:size(ripples.peaks,1)
                    ripple_ave(i,:) = eeg(round(ripples.peaks(i)*LfpSamplingrate)-Win:round(ripples.peaks(i)*LfpSamplingrate)+Win,:);
        end
        All_Ripple_Avg=double(mean(ripple_ave));
        Rippl_Matrix{1,shk}(CH,:)=All_Ripple_Avg;
    end
end

%%
figure('position',[200 115 1300 800])
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
        
        plot(Ripple_X,Ripple_Y,'color','k','linewidth',1);
        hold on
    end
end

set(gca,'YTick',[],'XTick',[]);set(gca,'color','w');
title([basename '  SWR profile']); 
saveas(gca,[basename '.SWRprofile.fig'])

end
