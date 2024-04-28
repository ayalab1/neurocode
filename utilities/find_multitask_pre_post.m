function pre_task_post = find_multitask_pre_post(env,varargin)
% Finds the indices of pre-task sleep/ task / post-task sleep from
% session.epochs. Input 'env' easily obtained by using load_epoch function
% to load session epochs. 

% input: 
%  env: cell array from epoch.name (epoch is dataframe from load_epoch).
% (optional)
% task_tag: cell array (default{'open_field','linear_track'}) enviornment 
%   designations that indicate your task. 
% post_sleep_flank: boolean (default: false), indicates post task sleep must follow task.
% pre_sleep_common: boolean (default: true), indicates if the first sleep
%   session should be used as reference for multitask sessions. 
% output: 
% pre_task_post: n x 3 matrix, where n is the number of tasks with
% pre/task/post.  
% 
% 
% LB 2024 ported from neuro_py find_multitask_pre_post

p = inputParser;
addParameter(p,'task_tag',{'open_field','linear_track'},@iscell);
addParameter(p,'post_sleep_flank', false,@isnumeric);
addParameter(p,'pre_sleep_common',true,@isbool);

parse(p,varargin{:});
task_tag = p.Results.task_tag;
post_sleep_flank = p.Results.post_sleep_flank;
pre_sleep_common = p.Results.pre_sleep_common;


    % Find the row indices that contain the search string in the specified column
    if isempty(task_tag)
        task_bool = ~contains(env, 'sleep', 'IgnoreCase', true);
    else
        task_bool = contains(env, task_tag, 'IgnoreCase', true);
    end
    sleep_bool = contains(env, 'sleep', 'IgnoreCase', true);

    task_idx = find(task_bool);
    task_idx(task_idx == 1) = [];
    sleep_idx = find(sleep_bool);

    pre_task_post = [];
    for task = task_idx
        temp = sleep_idx - task;
        pre_task = sleep_idx(temp < 0);
        post_task = sleep_idx(temp > 0);

        if isempty(post_task)
            warning(['no post_task sleep for task epoch ', num2str(task)]);
        elseif isempty(pre_task)
            warning(['no pre_task sleep for task epoch ', num2str(task)]);
        else
            pre_task_post = [pre_task_post; [pre_task(end), task, post_task(1)]];
        end
    end

    if isempty(pre_task_post)
        pre_task_post = [];
    end

    % search for epochs where the last epoch is 1 more than the first epoch
    if post_sleep_flank && ~isempty(pre_task_post)
        pre_task_post_ = [];
        for seq = pre_task_post'
            if seq(3) - seq(2) == 1
                pre_task_post_ = [pre_task_post_; seq'];
            end
        end
        pre_task_post = pre_task_post_;
    end

    % make the first pre task sleep the same pre task in subsequent tasks
    if pre_sleep_common && ~isempty(pre_task_post)
        pre_task_post_ = [];
        for seq = pre_task_post'
            pre_task_post_ = [pre_task_post_; [pre_task_post(1,1), seq(2), seq(3)]];
        end
        pre_task_post = pre_task_post_;
    end
        
end
