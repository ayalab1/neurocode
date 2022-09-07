% function [ratio_colab] = doulbe_labeling_ratio(file_a, file_b)
% This function finds the double-labeling of neurons under immunofluorescence
% by quantifying the ratio of colabeled neurons over the number of total
% number of cells.
% E.g. cell type a: tdTom; cell type b: brdu -> ratio_colab = ab / (a+b)

file_a = "C3-0.vsi_-_10x_05";
file_b = "C2-0.vsi_-_10x_05";

% set the scale for convertion: 
scale = 0.5725;

% set radius of detection of colabeled cell
r_detect = 15; %um
r_detect_pixel = r_detect*scale;

cells_coordinates_file_dir=pwd+"\";
% Extract coordinates
data_table_a = readtable(cells_coordinates_file_dir+file_a+".csv");
data_X_a = table2array(data_table_a(:,2));
data_Y_a = table2array(data_table_a(:,3));
data_table_b = readtable(cells_coordinates_file_dir+file_b+".csv");
data_X_b = table2array(data_table_b(:,2));
data_Y_b = table2array(data_table_b(:,3));
% number of cells
nCells_a=length(data_X_a);
nCells_b=length(data_X_b);

%% Plotting
figure

set(gca, 'YDir','reverse')

bw = imread(file_b+'.png');
imshow(bw)
hold on
scatter(data_X_a,data_Y_a,'b')
scatter(data_X_b,data_Y_b,'r*')
%%

% Find double-labeling cells
colabeled_index_a = [];
colabeled_index_b = [];
nColab = 0;
for i = 1:nCells_a
    for j = 1:nCells_b
        x1=data_X_a(i);y1=data_Y_a(i);
        x2=data_X_b(j);y2=data_Y_b(j);
        if (x1-x2)^2+(y1-y2)^2 < r_detect_pixel^2
            nColab = nColab+1;
            colabeled_index_a = [colabeled_index_a i];
            colabeled_index_b = [colabeled_index_b j];
            i
        end
    end
end

%% plot colabeled cells
figure

set(gca, 'YDir','reverse')

bw = imread(file_b+'.png');
imshow(bw)
hold on
scatter(data_X_a(colabeled_index_a),data_Y_a(colabeled_index_a),'b')
scatter(data_X_b(colabeled_index_b),data_Y_b(colabeled_index_b),'r*')
% end




