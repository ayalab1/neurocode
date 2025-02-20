clc, close all, clear all

dx = 0.8; % amount of blue and red at the beginning & end of the colormap
n = 32;  % number of entries in the rainbow colormap
[cmapRainbow, g, f]=rainbow_cmap(n, dx);
% % -----------------------------------------------
% % Display an image in bone and rainbow colormaps
load spine 
figure, image(X), colormap(bone), title('\bf  Image in Bone colormap');
figure, image(X), colormap(cmapRainbow), title('\bf Image in Rainbow colormap');

% % -----------------------------------------------
% % Plot transfer functions of rainbow colormap
figure, hold on
plot(g,cmapRainbow(:,1),'r','linewidth',2); %red
plot(g,cmapRainbow(:,2),'g','linewidth',2); %green
plot(g,cmapRainbow(:,3),'b','linewidth',2); %blue
% legend('R','G','B');
handle1=gca;
SetXLimYLimXTick(handle1,0,6,-0.1,1.1);
grid on
title('\bf Rainbow colormap transfer functions');
xlabel('\bf Scalar value g');
ylabel('\bf Transfer functions: R, G, B');


% % -----------------------------------------------
% % This program or any other program(s) supplied with it does not provide any warranty direct or implied.
% % This program is free to use/share for non-commercial purpose only. 
% % Kindly reference the work.
% % Author: Dr. Murtaza Khan
% % LinkedIn: http://www.linkedin.com/pub/dr-murtaza-khan/19/680/3b3
% % ResearchGate: http://www.researchgate.net/profile/Murtaza_Khan2/
% % -----------------------------------------------



