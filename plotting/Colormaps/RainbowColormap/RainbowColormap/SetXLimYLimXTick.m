% % function to set ticks of X-axis and Y-axis

function SetXLimYLimXTick(handle,xmin,xmax,ymin,ymax)

set(handle,'XLim',[xmin xmax]);
set(handle,'YLim',[ymin ymax]);

xtvec=[xmin ,get(gca,'XTick'), xmax];
ytvec=[ymin ,get(gca,'YTick'), ymax];

xtvec=unique(xtvec);
ytvec=unique(ytvec);

set(handle,'XTick',xtvec);
set(handle,'YTick',ytvec);

% % ---------------------------------------
% % This program or any other program(s) supplied with it does not provide any warranty direct or implied.
% % This program is free to use/share for non-commercial purpose only. 
% % Kindly reference the work.
% % Author: Dr. Murtaza Khan
% % LinkedIn: http://www.linkedin.com/pub/dr-murtaza-khan/19/680/3b3
% % ResearchGate: http://www.researchgate.net/profile/Murtaza_Khan2/
% % ---------------------------------------
