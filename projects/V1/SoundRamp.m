function [sound] = SoundRamp(startFrequency,stopFrequency,duration,samplingRate,gaussian)

% EXAMPLE:
% Play upward ramp:
% sound = SoundRamp(5000,20000,0.8,200000,1);
% handle = audioplayer(sound(:,2),1/mode(diff(sound(:,1)))); 
% playblocking(handle);

if ~exist('startFrequency','var'), startFrequency = 5000; end
if ~exist('stopFrequency','var'), stopFrequency = 10000; end
if ~exist('duration','var'), duration = 0.8; end
if ~exist('samplingRate','var'), samplingRate = 100000; end
if ~exist('gaussian','var'), gaussian = 1; end

if gaussian==0 % no gaussian fade-in or out
    dt = 1./samplingRate;
    t = (1:round(duration/dt))'*dt;
    f = linspace(startFrequency,stopFrequency,length(t))';
    y = cos(cumsum(f*dt*2*pi));
    sound = [t y];
elseif gaussian==1 % fade out only (avoid abrupt change in signal, resulting in a sound clicking effect)
    add = 0.1;
    width = add/5;
    stopFrequency1 = stopFrequency + add*((stopFrequency-startFrequency)/0.8);
    duration1 = duration + add;
    sound = SoundRamp(startFrequency,stopFrequency1,duration1,samplingRate,0);
    fade = double(InIntervals(sound(:,1),[0 duration]));
    in = sound(:,1)>duration;
    fade(in) = gaussmf(sound(in,1),[width duration]);
    sound(:,2) = sound(:,2).*fade;
elseif gaussian==2 % fade in and fade out
    add = 0.1;
    width = add/5;
    stopFrequency1 = stopFrequency + add*(((stopFrequency-startFrequency))/0.8);
    startFrequency1 = startFrequency + add*(((startFrequency-stopFrequency))/0.8);
    duration1 = duration + add*2;
    sound = SoundRamp(startFrequency1,stopFrequency1,duration1,samplingRate,0);
    fade = double(InIntervals(sound(:,1),[0 duration]+add));
    in = sound(:,1)<add;
    fade(in) = gaussmf(sound(in,1),[width add]);
    in = sound(:,1)>duration+add;
    fade(in) = gaussmf(sound(in,1),[width duration+add]);
    sound(:,2) = sound(:,2).*fade;
end

