function pixel_distance = maze_distance_gui(video_path)
% maze_distance_gui: manually get distance between points in pixels
%
% This fuction was made in order to get the ratio to convert tracking
% points from pixels to cm.
%
% Input: 
%   video_path: full path to video
%
% Output:
%   pixel_distance: distance between clicked points in pixels
%
% Ryan H

% check if video exists
if ~exist(video_path,'file')
   error([video_path,'  video does not exist']) 
end

% read video
vid_obj = VideoReader(video_path);

% read first 10 seconds of video frames
frames = read(vid_obj, [1, round(vid_obj.FrameRate*10)]);
% init matrix to store flattened frames
grey_frames = zeros(vid_obj.Height, vid_obj.Width, size(frames,4));
for i = 1:size(frames, 4)
    grey_frames(:,:,i) = rgb2gray(frames(:,:,:,i));
end
% take averge of 10 seconds
grey_frames_avg = mean(grey_frames, 3);

% plot average frame
fig = figure;
imagesc(grey_frames_avg)
hold on;
colormap('gray');
axis('image')
title('click on 2 key points and hit "enter" ')

% let the user click around the coordinates
corners = [];
while true
    [X,Y]=ginput(1);
    if isempty(X)
        break
    end
    corners = [corners;[X,Y]];
    plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
    i=i+1;
end
close(fig)

% calculate pixel distance
pixel_distance = pdist(corners,'euclidean');
end

