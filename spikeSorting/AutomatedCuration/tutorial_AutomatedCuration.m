% A little explanation for how to get the automated curation up and running!

%% ============================= INSTALLATION =============================
% First, if you are doing this for the very first time, you need to make a 
% python environment with the required components installed. Here is how to
% do this on Windows:
%
% Run (in a cosole or anaconda prompt): [NOT matlab]
% conda create --name AutomatedCurator
% conda activate AutomatedCurator  
% conda install -c fastchan fastai anaconda      % agree to install when prompted
% pip install phylib      

%% =============================== PREP DATA ===============================

% If you want to cluster-cut on your local machine's SSD drive, now is the 
% time to copy your phy files and the .dat file to your local drive, as you
% normally, would before this step.

% Open phy and use the plugin: "Recluster Local PCAs all"
% This will take some time, but this is a very important step of 
% oversplitting to make sure clusters are sufficiently clean to be 
% run through the automated curator.


% You absolutely need your .xml file in the same folder.
% Note - Launch this within Kilosort folder 
clustering_path = pwd; 
[parentFolder] = fileparts(clustering_path);
[~,basename] = fileparts(parentFolder);
neurosuite_path = parentFolder; % the files will be saved with the dat file.
datFile = fullfile(parentFolder,[basename '.dat']);
% The reason for saving the files with the dat file is that the 
% Automated Curator expects the filenames to be the same name as the name of the folder (and the .xml)

% Launch the convertion of the phy format into the Neurosuite format that the 
% Automated Curator needs
% NOTE - Launch this within Kilosort folder 
Phy2Neurosuite(clustering_path,parentFolder,datFile);

%% ======================= LAUNCH AUTOMATED CURATION =======================

% From within your conda environment (conda activate AutomatedCurator),
% launch (for each shank), the following:

% Prepare a script to do an automated call:
channelShanks = double(readNPY(fullfile(clustering_path, 'channel_shanks.npy')));

github_path = 'C:\Users\Cornell\Documents\GitHub\';
if ~exist(github_path,'dir')
    github_path = 'D:\github\';
end
fid = fopen(fullfile(neurosuite_path,'launch_curation.bat'),'w');
fwrite(fid, sprintf('%s\n', '@echo_on'));
fwrite(fid, sprintf('%s\n', 'call conda deactivate'));
fwrite(fid, sprintf('%s\n', 'call conda activate AutomatedCurator'));
lineStart = ['python ' fullfile(github_path,'neurocode','spikeSorting','AutomatedCuration','Automated-curation','running_AI_pipeline_bash.py') ' ' strrep(fullfile(neurosuite_path,basename),'\','\\')];

cluster_info_table = importdata(fullfile(clustering_path,'cluster_info.tsv')); 
shankID = cluster_info_table.data(:,end); 
nShanks = max(shankID);
for i=unique(shankID(:))'
    line = [lineStart ' ' num2str(i) ' ' num2str(sum(channelShanks==i))];
    fwrite(fid, sprintf('%s\n\n', line));
end
fclose(fid);

% Run the file! 
%NOTE - run outside of kilosort folder but inside day folder (where .res, .fet,.spk files are)
evalc(['!' fullfile(neurosuite_path,'launch_curation.bat')]); % Here, we are trying to run the file from matlab. 
% If you cannot run the file from matlab, run the newly created "launch_curation.bat" file manually. You may need to click "Enter" to help it along...

%% ======================= EXPORT BACK TO PHY =======================



%NOTE - Run inside kilosort folder (ensure clustering_path, neurosuite_path is run
%INSIDE kilosort folder, if not then below will look for .clu in wrong dir)
UpdatePhyFromNeurosuite(clustering_path,neurosuite_path);
%if you have issues finding cluster_info.tsv - open up kilosort folder in
%phy and save once and cluster_info.tsv will be created

% Enjoy!


