function [intstimestamps_samples,intstimestamps] = make_pre_task_post(basePath,session_sequence)
% generate pre/task/post structure

% session_sequence: type of sub-sessions in order of recording. 1=
% preSleep, 2= task, 3=postSleep

% Farnaz

% Example :
% basePath ='F:\OML\OML18\day1';
% session_sequence=[1 2 2 3];

warning('this function does not work well with discontinuos task epochs')

[~,recordingname] = fileparts(basePath);
load([recordingname '.MergePoints.events.mat']);
pre_idx_s=[];task_idx_s=[];post_idx_s=[];
pre_idx=[];task_idx=[];post_idx=[];

if nargin > 1
    if length(session_sequence)~=size(MergePoints.foldernames,2)
        error(['Number of the folders must be equal to number of the folders in Mergepoints = ',...
            [num2str(size(MergePoints.foldernames,2))], ' (MergePoints.foldernames)']);
    end
    pre=find(session_sequence==1);
    task=find(session_sequence==2);
    post=find(session_sequence==3);
    
    pre_idx_s=MergePoints.timestamps_samples(pre,:);
    task_idx_s=MergePoints.timestamps_samples(task,:);
    post_idx_s=MergePoints.timestamps_samples(post,:);
    
    pre_idx=MergePoints.timestamps(pre,:);
    task_idx=MergePoints.timestamps(task,:);
    post_idx=MergePoints.timestamps(post,:);
else
    for cond=1:size(MergePoints.foldernames,2)
        
        FN=MergePoints.foldernames{1,cond};
        pre=strfind(lower(FN),'pre');
        post=strfind(lower(FN),'post');
        
        if isempty(pre)==0
            pre_idx_s=[pre_idx_s;MergePoints.timestamps_samples(cond,:)];
            pre_idx=[pre_idx;MergePoints.timestamps(cond,:)];
        elseif isempty(post)==0
            post_idx_s=[post_idx_s;MergePoints.timestamps_samples(cond,:)];
            post_idx=[post_idx;MergePoints.timestamps(cond,:)];
        else
            task_idx_s=[task_idx_s;MergePoints.timestamps_samples(cond,:)];
            task_idx=[task_idx;MergePoints.timestamps(cond,:)];
        end
    end
end

intstimestamps_samples.pre_NonMerged=pre_idx_s;
intstimestamps_samples.task_NonMerged=task_idx_s;
intstimestamps_samples.post_NonMerged=post_idx_s;

intstimestamps.pre_NonMerged=pre_idx;
intstimestamps.task_NonMerged=task_idx;
intstimestamps.post_NonMerged=post_idx;

intstimestamps_samples.pre=[pre_idx_s(1,1) pre_idx_s(end,2)];
intstimestamps_samples.task=[task_idx_s(1,1) task_idx_s(end,2)];
intstimestamps_samples.post=[post_idx_s(1,1) post_idx_s(end,2)];

intstimestamps.pre=[pre_idx(1,1) pre_idx(end,2)];
intstimestamps.task=[task_idx(1,1) task_idx(end,2)];
intstimestamps.post=[post_idx(1,1) post_idx(end,2)];

behavEpochs.int_samples=intstimestamps_samples;
behavEpochs.int=intstimestamps;
save([recordingname '.behavEpochs.mat'],'behavEpochs')
end