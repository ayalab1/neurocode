% Read the fluorescence level of all cells in all the images in the folder
% This version uses the surrounding pixels in a rectangle to define cell
% fluorescence
%% Setting
roi_radius = 3;

%%
error("Please switch the pwd to the png folder.")
%% Read image (red channel)
imagefiles = dir('*.png');
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
    currentfilename = imagefiles(ii).name;
    currentimage = imread(currentfilename);
    [folder, baseFileNameNoExt{ii}, extension] = fileparts(currentfilename);
    images{ii} = currentimage;
end


%%
% change to the PNG folder
for ii=1:nfiles
    % Read in original RGB image.
    rgbImage = imread(imagefiles(ii).name);
    % Extract color channels.
    redChannel = rgbImage(:,:,1); % Red channel
    greenChannel = rgbImage(:,:,2); % Green channel
    blueChannel = rgbImage(:,:,3); % Blue channel
    % Create an all black channel.
    allBlack = zeros(size(rgbImage, 1), size(rgbImage, 2), 'uint8');
    % Create color versions of the individual color channels.
    just_red = cat(3, redChannel, allBlack, allBlack);
    just_green = cat(3, allBlack, greenChannel, allBlack);
    just_blue = cat(3, allBlack, allBlack, blueChannel);
    
    % Save red channels
    images_red{ii} = redChannel;
    % Plot separate channels
    % % Recombine the individual color channels to create the original RGB image again.
    % recombinedRGBImage = cat(3, redChannel, greenChannel, blueChannel);
    % % Display them all.
    % subplot(3, 3, 2);
    % imshow(rgbImage);
    % fontSize = 20;
    % title('Original RGB Image', 'FontSize', fontSize)
    % subplot(3, 3, 4);
    % imshow(just_red);
    % title('Red Channel in Red', 'FontSize', fontSize)
    % subplot(3, 3, 5);
    % imshow(just_green)
    % title('Green Channel in Green', 'FontSize', fontSize)
    % subplot(3, 3, 6);
    % imshow(just_blue);
    % title('Blue Channel in Blue', 'FontSize', fontSize)
    % subplot(3, 3, 8);
    % imshow(recombinedRGBImage);
    % title('Recombined to Form Original RGB Image Again', 'FontSize', fontSize)
    % % Set up figure properties:
    % % Enlarge figure to full screen.
    % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    % % Get rid of tool bar and pulldown menus that are along top of figure.
    % % set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    % % Give a name to the title bar.
    % set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
end

%%
error("Please switch the pwd to the processed data folder.")
%% Extract cell fluorescence intensity
for ii=1%:nfiles
    data_table = readtable(strcat(baseFileNameNoExt{ii},'.csv'));
    data_X = table2array(data_table(:,2));
    data_Y = table2array(data_table(:,3));
    nCells = length(data_X);
    size_image = size(images_red{ii});
    images_temp = images_red{ii};
    
    for i_cell=1:nCells
        left_bound = max(0,round(data_X(i_cell)-roi_radius));
        right_bound = min(size_image(2),round(data_X(i_cell)+roi_radius));
        up_bound = min(size_image(1),round(data_Y(i_cell)+roi_radius));
        down_bound = max(0,round(data_Y(i_cell)-roi_radius));
        intensity_list{ii}{i_cell} = sum(sum(images_red{ii}(down_bound:up_bound,left_bound:right_bound)));
        images_temp(down_bound:up_bound,left_bound:right_bound) = intensity_list{ii}{i_cell}/50;
    end
    
    %Plot rois
    figure
    hold on
    image(images_red{ii})
    
    for i_cell=1:nCells
        left_bound = max(0,round(data_X(i_cell)-roi_radius));
        right_bound = min(size_image(2),round(data_X(i_cell)+roi_radius));
        up_bound = min(size_image(1),round(data_Y(i_cell)+roi_radius));
        down_bound = max(0,round(data_Y(i_cell)-roi_radius));
        intensity_list{ii}{i_cell} = sum(sum(images_red{ii}(down_bound:up_bound,left_bound:right_bound)));
        rectangle('Position',[left_bound down_bound 6 6])
    end
    text(data_X,data_Y,intensity_list{ii})
end

%% Concat all the cell intensity
intensity_list_all_cells = [];
for ii=1:nfiles
    intensity_list_all_cells = [intensity_list_all_cells, cell2mat(intensity_list{ii})];
end
cd ..
[~,brain_name,~]=fileparts(pwd);
figure
t = title(strcat("Distribution of fluorescence intensity of ",brain_name))
set(t, 'Interpreter', 'none')
ylabel("Count")
xlabel("Fluorescence intensity")
hold on
h = histogram(intensity_list_all_cells,50);
h.EdgeColor = "None"
saveas(gcf,strcat("Distribution of fluorescence intensity of ",brain_name,'.fig'))