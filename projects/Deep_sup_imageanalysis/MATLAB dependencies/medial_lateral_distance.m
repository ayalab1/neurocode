function [medial_lateral_position, normalized_medial_lateral_position] = medial_lateral_distance(file)
% This function calculate the distance of the projection of each cell onto
% the superficial border of CA1 pyramidal layer to the most medial point of
% CA1, and then normalize these distances.

% Instruction: please run cell_depth.m using main.m first to get all the
% necessary data. Then, find the most medial and lateral point of CA1 in
% ImageJ.

%   Input
%   file                                name of the file to be processed

%   Output
%   medial_lateral_position             distance along the curve (the superficial border of CA1 pyramidal layer) from each
%                                           projection foot of each cell to the most medial point
%   normalized_medial_lateral_position  normalized
%                                           (medial_lateral_position) s.t. the most lateral point has distance "1"
%% Load
load(file+".mat")
%% Find the projection of the most medial and lateral point onto the superficial border
medialFileID = fopen(file+"_medial.txt",'r');
laterlFileID = fopen(file+"_lateral.txt",'r');
formatSpec = "%f";
medial_point= fscanf(medialFileID,formatSpec)';
lateral_point = fscanf(laterlFileID,formatSpec)';
[medial_pFoot,~,medial_position] = distance2curve([line_X_down,line_Y_down],[[medial_point(1)],[medial_point(2)]]);
[lateral_pFoot,~,lateral_position] = distance2curve([line_X_down,line_Y_down],[[lateral_point(1)],[lateral_point(2)]]);

%% Find the location of each point (each cell projection and the medial, lateral point projection)
% define a vector of position of cells, normalized by the length of the
% border
position_cell = zeros(1,nCells);
for i = 1:nCells
    [~,~,position_cell(i)] = distance2curve([line_X_down,line_Y_down],[data_X(i),data_Y(i)],'spline');
end

% Re-normalize the position of the cells based on the medial and lateral
% point projection position
if medial_position > lateral_position % starts with the lateral part, inversed!
    normalized_medial_lateral_position = 1 - (position_cell - lateral_position)/(medial_position - lateral_position);
else
    normalized_medial_lateral_position = (position_cell - medial_position)/(lateral_position - medial_position);
end

%%
% calculate the real position (distance along the border to the medial
% point) by using the scaling factor
length_border = arclength(line_X_down,line_Y_down,'spline')/scale;
if medial_position > lateral_position % starts with the lateral part, inversed!
    medial_lateral_position = (1 - position_cell) * length_border;
else
    medial_lateral_position = position_cell * length_border;
end

%% plot
figure
hold on
axis equal
plot(line_X_down,line_Y_down);
plot(line_X_up,line_Y_up);
scatter(data_X,data_Y,[],[0.3,0.5,0.8])
scatter([pHead(:,1);pFoot(:,1)],[pHead(:,2);pFoot(:,2)],'.','black')
scatter(medial_pFoot(1),medial_pFoot(2),"rx")
scatter(lateral_pFoot(1),lateral_pFoot(2),"rx")
legend('SP-SR','SP-SO','cell','projection','medial and lateral end')

ax = gca;
ax.YDir = 'reverse';
dim = [.2 .6 .3 .3];
str = 'deep';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
dim = [.2 .01 .3 .3];
str = 'superficial';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
dim = [.5 .6 .3 .3];
annotation('textbox',dim,'String',file,'FitBoxToText','on');
xLabel = strcat('Unit of plot: pixel; Scale of original figure: ',string(scale),' _ pixle/um');
xlabel(xLabel)
xlim([200 1500])
ylim([0 700])

%% save file
save(file);
disp(file)
disp('done')
end
