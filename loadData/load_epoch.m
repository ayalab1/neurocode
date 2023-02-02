function epochs = load_epoch(varargin)
% load_epoch: loads nice table of session sub-epochs from basename.session
%
% Input:
%   basepath - path to session folder
%
% Output:
%   epochs - table with session epoch metadata
%
% Ryan Harvey 2023

p=inputParser;
addParameter(p,'basepath',pwd);
parse(p,varargin{:});
basepath_ = p.Results.basepath;
basename = basenameFromBasepath(basepath_);

% load session file
load(fullfile(basepath_,[basename,'.session.mat']),'session')

% iter over epochs and pull out metadata
for ep = 1:length(session.epochs)
    
    name{ep,1} = session.epochs{ep}.name;
    startTime(ep,1) = session.epochs{ep}.startTime;
    stopTime(ep,1) = session.epochs{ep}.stopTime;
    
    if isfield(session.epochs{ep},'environment')
        environment{ep,1} = session.epochs{ep}.environment;
    else
        environment{ep,1} = [];
    end
    
    if isfield(session.epochs{ep},'behavioralParadigm')
        behavioralParadigm{ep,1} = session.epochs{ep}.behavioralParadigm;
    else
        behavioralParadigm{ep,1} = [];
    end
    if isfield(session.epochs{ep},'manipulation')
        manipulation{ep,1} = session.epochs{ep}.manipulation;
    else
        manipulation{ep,1} = {};
    end
    
    if isfield(session.epochs{ep},'stimuli')
        stimuli{ep,1} = session.epochs{ep}.stimuli;
    else
        stimuli{ep,1} = [];
    end
    
    if isfield(session.epochs{ep},'notes')
        notes{ep,1} = session.epochs{ep}.notes;
    else
        notes{ep,1} = [];
    end
    
    basepath{ep,1} = basepath_;
end

% place into table
epochs = table(name,startTime,stopTime,environment,behavioralParadigm,...
    manipulation,stimuli,notes,basepath);

end