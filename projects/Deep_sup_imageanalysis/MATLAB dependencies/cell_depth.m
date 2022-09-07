function [pFoot,sign,depth_abs,data_X,data_Y,t_a,line_X_down,line_Y_down,line_X_up,line_Y_up,file,scale,SP_thickness,depth_norm] = cell_depth(file)
% version2 : 20210224 update: calculate the normalized depth
% 20210214 save 2 lines coordinates but not calculate ration

%This function finds the projection from the cell to the SP-SR border, and
%then find the intersection of the projection line and the SP-SO border.
%Then the SP thickness is defined by the Euclidean distance between pFoot
%pHead. The absolute depth is defined by the distance between the cell and pFoot.
%Then, the normalized depth is defined by (absolute depth) / (SP thickness)

%Instruction: please process the data in ImageJ using the macros. Copy the
%directory of the saved cvs file, txt file as the input here.
%   pFoot: coordinates of projection to the SP-SR border
%   pHead: coordinates of the intersection of the projection and the SP-SO border
%   depth_abs: absolute distance to border (plus->deep)
%   depth_norm: normalized depth
%   SP_thickness: thickness of the pyramidal cell layer
%   data_: coordinates of the data points
%   line_: coordinates of SP-SO (up)/ SP-SR(down) border
%   t_a: fractional arc length along the interpolating curve to that point.

%setting the scale in OLYvis to the max of 200um. (This scale only work on Yunchang's PC)
%the scales for each sections are stored in temp_data

scale = 0.5725;
% read cells coordinates
% cells_coordinates_file_dir = "D:\微云同步盘\George's Cloud\Abroad\NYU\工作\Buzsaki Lab\Projects\Roman\Development\Histology analysis\ImageJ results\";
cells_coordinates_file_dir=pwd+"\";
data_table = readtable(cells_coordinates_file_dir+file+".csv");
data_X = table2array(data_table(:,2));
data_Y = table2array(data_table(:,3));
nCells=length(data_X);
%% load borders coordinate

downfileID = fopen(file+"_down.txt",'r');
upfileID = fopen(file+"_up.txt",'r');
upfileID
formatSpec = "%f";
downCoordinates= fscanf(downfileID,formatSpec);
line_X_down = downCoordinates(1:length(downCoordinates)/2);
line_Y_down = downCoordinates(length(downCoordinates)/2+1:end);

upCoordinates = fscanf(upfileID,formatSpec);
line_X_up = upCoordinates(1:length(upCoordinates)/2);
line_Y_up = upCoordinates(length(upCoordinates)/2+1:end);

%% Smooth the boarders with splines
% smoothing the down boarder
x = line_X_down';
y = line_Y_down';
% Cubic spline data interpolation
t = 1:numel(x);
xy = [x;y];
pp = spline(t,xy);
tInterp = linspace(1,numel(x));
xyInterp = ppval(pp, tInterp);
line_X_down=xyInterp(1,:)';
line_Y_down=xyInterp(2,:)';

% smoothing the up boarder
x = line_X_up';
y = line_Y_up';
% Cubic spline data interpolation
t = 1:numel(x);
xy = [x;y];
pp = spline(t,xy);
tInterp = linspace(1,numel(x));
xyInterp = ppval(pp, tInterp);
line_X_up=xyInterp(1,:)';
line_Y_up=xyInterp(2,:)';

%% Find pFoot, pHead, abs depth, and SP_thickness
[pFoot,depth_abs,t_a] = distance2curve([line_X_down,line_Y_down],[data_X,data_Y]);
[sign] = aboveOrBelow(pFoot,data_X,data_Y);
% extend the line between the cell and pFoot so that it intersects with the
% upper border (SP-SO)
%%

pHead=NaN(nCells,2);%rows: coordinates; collumns: x,y

for i=1:nCells
    intersect_boo = 0;%boolean to indicate when the intersection is first found
    x1=pFoot(i,1);x2=data_X(i);y1=pFoot(i,2);y2=data_Y(i);
    while intersect_boo == 0
        [x1,x2,y1,y2] = expand_line(x1,x2,y1,y2);
        [p] = InterX([x1,x2;y1,y2],[line_X_up';line_Y_up']);
        if(x1==x2&&y1==y2) %cell is on the down border, then pHead is the Euclidean distance to up line
            [p,~,~] = distance2curve([line_X_up,line_Y_up],[x1,y1]);
            disp(strcat('in file: ',file))
            disp('cell on border')
            pHead(i,1)=p(1);
            pHead(i,2)=p(2);
            intersect_boo=1;
        elseif abs(x1)+abs(x2)+abs(y1)+abs(y2)>10000
            % kill the program is not working
            disp(file)
            %plot
            figure
            hold on
            xlim([0 1000])
            ylim([0 1000])
            axis equal
            title(strcat('Not intersecting: ',file))
            plot(line_X_down,line_Y_down);
            plot(line_X_up,line_Y_up);
            scatter(data_X,data_Y,[],[0.3,0.5,0.8])
            scatter(pFoot(:,1),pFoot(:,2),'.','black')
            legend('SP-SR','SP-SO','cell','projection')
            ax = gca;
            ax.YDir = 'reverse';
            dim = [.2 .6 .3 .3];
            str = 'deep';
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            dim = [.2 .01 .3 .3];
            str = 'superficial';
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            xLabel = strcat('Unit of plot: pixel; Scale of original figure: ',string(scale),' _ pixle/um');
            xlabel(xLabel)
            %error message
            error('not intersecting')
        elseif ~isempty(p)
            pHead(i,1)=p(1,1);
            pHead(i,2)=p(2,1);
            intersect_boo=1;
        end
    end
end

%% Calculate final results, plot, and save
SP_thickness = NaN(nCells,1);
for i =1:nCells
    x1 = pFoot(i,1);x2 = pHead(i,1);y1 = pFoot(i,2);y2 = pHead(i,2);
    SP_thickness(i)=sqrt((x1-x2)^2+(y1-y2)^2);
end
depth_abs = -sign.*depth_abs;
depth_norm = depth_abs./SP_thickness;
%set everything to the scale of the sample. positive ->deep
depth_abs = depth_abs./scale;
SP_thickness=SP_thickness./scale;

%% plot
figure
hold on
axis equal
plot(line_X_down,line_Y_down);
plot(line_X_up,line_Y_up);
scatter(data_X,data_Y,[],[0.3,0.5,0.8])
scatter(pHead(:,1),pHead(:,2),'.','black')
scatter(pFoot(:,1),pFoot(:,2),'.','black')
legend('SP-SR','SP-SO','cell','projection')
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
%save file
save(file);
disp(file)
disp('done')

end
