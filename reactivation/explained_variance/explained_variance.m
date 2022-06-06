function [EV, rEV] = explained_variance(spikes,intervals,varargin)
% Calculate explained variance of the firing rate correlation matrix of
% several behavioral epochs
% 
% INPUTS
% spikes    struct  spikes.cellinfo struct
% intervals	cell    1x3 array with intervals (one or multiple)to be
%                   used as pre, beh and post intervals
% 
% For details see Kudrimoti, Barnes and McNaughton 1999
% Adapted to neurocode by aza, 2021

p = inputParser;

addParameter(p,'UIDs',[],@isnumeric);
addParameter(p,'bin',.05,@isnumeric); % Time bin for firing rate vectors (sec)
addParameter(p,'win',0,@isnumeric); % Time window for padding around Intervals
parse(p,varargin{:})

UIDs = p.Results.UIDs;
bin = p.Results.bin;
win = p.Results.win;

%test input
if isempty(spikes)
    error(['No spikes data']);
end

if isempty(UIDs)
    UIDs = spikes.UID;
end

%adjust intervals with provided padding window
for i = 1:size(intervals,1)
    intervals{i}(1,:) = intervals{i}(1,:) - win;
    intervals{i}(2,:) = intervals{i}(1,:) + win;
end

%% Calculate firing rate vectors pre/post
incl = true(1,length(UIDs));
for ep = 1:size(intervals,1)
    ratemaps{ep} = [];
    for i = 1:size(intervals{ep},1)
        for u = 1:length(UIDs)
            tmp(:,u) = histcounts(spikes.times{UIDs(u)},...
                intervals{ep}(i,1):bin:intervals{ep}(i,2));
            
%             % Exclude cells that fired zero spikes
%             if sum(theserates(:,u)) < 1
%                 incl(u) = false;
%             end
        end
        ratemaps{ep} = cat(1,ratemaps{ep},tmp);
        clear tmp
    end
    maplen(ep) = size(ratemaps{ep},1);
end
minr = min(maplen);
prelen = size(ratemaps{1},1);

% Exclude cells with zero spikes in relevant interval
incl = incl & sum(ratemaps{1}(prelen-minr+1:prelen,:),1)>0;
incl = incl & sum(ratemaps{2}>0);
incl = incl & sum(ratemaps{3}(1:minr,:),1)>0;

%% Calculate correlation matrices, forward-, and reverse explained variance
Rpre = cosmo_corr(ratemaps{1}(prelen-minr+1:prelen,incl)); 
RBeh = cosmo_corr(ratemaps{2}(:,incl));
Rpost = cosmo_corr(ratemaps{3}(1:minr,incl));

r_BehPost = corrcoef(RBeh(triu(RBeh,+1)~=0),Rpost(triu(Rpost,+1)~=0)); r_BehPost = r_BehPost(1,2);
r_BehPre = corrcoef(Rpre(triu(Rpre,+1)~=0),RBeh(triu(RBeh,+1)~=0)); r_BehPre = r_BehPre(1,2);
r_PrePost = corrcoef(Rpre(triu(Rpre,+1)~=0),Rpost(triu(Rpost,+1)~=0)); r_PrePost = r_PrePost(1,2);

EV = ((r_BehPost - r_BehPre*r_PrePost) / sqrt((1-r_BehPre^2)*(1-r_PrePost^2)))^2;

rEV = ((r_BehPre - r_BehPost*r_PrePost) / sqrt((1-r_BehPost^2)*(1-r_PrePost^2)))^2;

end