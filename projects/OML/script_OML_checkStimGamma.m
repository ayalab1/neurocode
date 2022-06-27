% function [z,basepath] = script_OML_checkStimGamma(basepath,varargin)
basepath = pwd;
cd(basepath);
figure
pulses = getStruct(basepath,'pulses');
session = getStruct(basepath,'session');
nChannels = session.extracellular.nChannels;
p = pulses.timestamps(:,1);
d = diff([0;p]);
q = p(d>1);

for i=1:length(q)
    interval = [-2 5]+q(i,1);
    lfp = GetAyaLFP(0:nChannels-1,'restrict',interval);
    for j=1:nChannels
        filtered = FilterLFP(lfp(:,[1 j+1]),'passband',[80 120]);
        [~,a] = Phase(filtered);
        mCell{i}(j,:) = a(:,2);
    end
    display([num2str(i) '/' num2str(length(q))]);
end
%%

lengths = cellfun(@(x) size(x,2),mCell);
mCell(lengths~=max(lengths)) = [];

leeway = 50;
x = linspace(interval(1)-q(i),interval(2)-q(i),size(a,1));
mm = cat(3,mCell{:});
mm = mm(:,leeway:end-leeway,:); 
x = x(leeway:end-leeway);
m = nanmedian(mm,3);
z = zscore(m,[],2);

clf

PlotColorMap(Smooth(z,[0,1250*0.2]),'x',x);
set(get(colorbar,'YLabel'),'String','gamma power (z-units)');
ylabel('channel ID');
xlabel('time from stimulation onset (s)');
set(gca,'box','off','fontsize',15);
clim([-1 1]*1.5);
title(basepath)
drawnow
saveas(gcf,fullfile(basepath,'StimGammaPowerAllChannels.png'));
saveas(gcf,fullfile(basepath,'StimGammaPowerAllChannels.fig'));
% clf
% semplot(x,z);