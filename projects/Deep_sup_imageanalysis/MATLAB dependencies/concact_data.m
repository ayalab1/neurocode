%% Concatinate all the cell depths within a folder
% This will load all data from the current directory. Be careful!

inpath = pwd;
files = dir(fullfile(strcat(inpath,'/*.mat')));
nFiles = length(files);%
depth_total=[]; % uncomment if var doesn't exist
for i=1:nFiles
    load(files(i).name)
    depth_total = cat(1,depth_total,depth_norm);
end
%% concat specific data in the workspace
depth_total_cat = cat(1,depth_total,depth);

%% Concat all the medial lateral positions together (normalized)
inpath = pwd;
files = dir(fullfile(strcat(inpath,'/*.mat')));
nFiles = length(files);%
medial_lateral_position_total=[]; % uncomment if var doesn't exist
for i=1:nFiles
    load(files(i).name)
    medial_lateral_position_total = cat(1,medial_lateral_position_total,normalized_medial_lateral_position');
end
%% plot
figure
histogram(medial_lateral_position_total,15)
xlabel("Medial - lateral position (Normalized)")
ylabel("Count")
t = title("E14_2M3 Medial - lateral position distribution")
set(t, 'Interpreter', 'none')

%% concat specific data in the workspace
medial_lateral_position_total_cat = cat(1,medial_lateral_position_total_cat,medial_lateral_position_total);
%% plot
figure
histogram(medial_lateral_position_total_cat,15)
xlabel("Medial - lateral position (Normalized)")
ylabel("Count")
t = title("E13 Medial - lateral position distribution")
set(t, 'Interpreter', 'none')

