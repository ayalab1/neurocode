function good_idx = manual_trackerjumps(ts,x,y,StartofRec,EndofRec,basepath,varargin)
% manual_trackerjumps: Allows you to manually cut out xy coordinates that
% are outside the bounds of your maze. These can be caused by unplugs or
% if the rat jumps out. If you do not remove these points, your ratemap
% will be messed up. 
%
% Input:
%       ts: timestamps 
%       x: x coord
%       y: y coords
%       StartofRec: timestamps indicating the start points of your event
%       EndofRec: timestamps indicating the end points of your event
%       darkmode: makes figure dark mode (true)
%       restic_dir: restrict to within shape (1) or to outside shape (0)
%       axis_equal: option to make x and y axis on the same scale
%       alpha: transparancy option of points (.2)
%       add_scatter: option to include scattered points (true)
%
% dependencies: darkBackground, basenameFromBasepath
%
% Ryan E Harvey (2018)


p = inputParser;
addParameter(p,'darkmode',true,@islogical);
addParameter(p,'restic_dir',1,@isnumeric);
addParameter(p,'axis_equal',true,@islogical);
addParameter(p,'alpha',.2,@isnumeric);
addParameter(p,'add_scatter',true,@islogical);

parse(p,varargin{:});
darkmode = p.Results.darkmode;
restic_dir = p.Results.restic_dir;
axis_equal = p.Results.axis_equal;
alpha = p.Results.alpha;
add_scatter = p.Results.add_scatter;

savets=[];
for i=1:length(StartofRec)
    % index out each event
    xtemp=x(ts>=StartofRec(i) & ts<=EndofRec(i));
    ytemp=y(ts>=StartofRec(i) & ts<=EndofRec(i));
    tstemp=ts(ts>=StartofRec(i) & ts<=EndofRec(i));
    % use the gui to cut out points
    [~,~,in]=restrictMovement(xtemp,ytemp,restic_dir,darkmode,axis_equal,alpha,add_scatter);
    % save the ts where the tracker error exists
    savets=[savets,tstemp(in)];
end
% locate the index for each tracker error
good_idx=ismember(ts,savets);
% save that index to your session folder so you won't have to do this again
% each time you run your data
basename = basenameFromBasepath(basepath);
save(fullfile(basepath,[basename,'.restrictxy.mat']),'good_idx')
end

function [x,y,in]=restrictMovement(x,y,direction,darkmode,axis_equal,alpha,add_scatter)
% restrictMovement allows you to draw a line around xy coordinates in order
% to eliminate certain points you don't want...ie when the rat jumps out of
% maze or tracker errors. 
%
% Also, I have included a direction input argument so you have either
% restrict outside or inside points. This is valuable if you have a maze
% like a circular track where the rat could jump out away or towards the
% center of the maze. 
%
% Input         x,y: coordinates 
%         direction: 1 (default) to remove outside points; 0 to remove inside points
%         
%
% Output        x,y: retained coordinates
%                in: logical of which coordinates were retained (so you can index ts etc.)
%       
%
% Ryan Harvey

% check inputs
if nargin<3
    direction=1;
end

% set up figure
if darkmode
    fig=figure;plot(x,y,'Color',[1,1,1,alpha]);hold on
    if add_scatter
        scatter(x,y,3,'w','filled');hold on
    end
    title('Click around the points you want to keep')
    xlabel('X')
    ylabel('Y')
    axis tight
    if axis_equal
        axis equal
    end
    darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
else
    fig=figure;plot(x,y,'Color',[0,0,0,alpha]);hold on
    if add_scatter
        scatter(x,y,3,'k','filled');hold on
    end
    title('Click around the points you want to keep')
    xlabel('X')
    ylabel('Y')
    axis tight
    if axis_equal
        axis equal
    end
end
disp('PRESS "ENTER" TO EXIT')
i=1;

% let the user click around the coordinates
while true
    [X,Y]=ginput(1);
    if isempty(X)
        break
    end
    corners(i,:)=[X,Y];
    plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
    i=i+1;
end

% remove points outside or inside the shape
in=inpolygon(x,y,corners(:,1),corners(:,2));
if direction==1
    x=x(in);
    y=y(in);
else
    x=x(~in);
    y=y(~in);
    in=~in;
end
close
end