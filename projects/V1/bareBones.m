
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movieDurationSecs = 2; % how long the display shows each stimulus
stimuli = [0 22.5 45  67.5  90  112.5  135 157.5]; %stimulus angles: 8 stimuli; -1 is black
nStimuli = length(stimuli);
nRepetitions = 120;
gratingsize = 800;
% res is the total size of the patch in x- and y- direction, i.e., the
% width and height of the mathematical support:
res = [gratingsize gratingsize];


% pre-set the pseudorandom stimulus order:
angles = [];
for j=1:nRepetitions
    theseAngles = stimuli(randperm(numel(stimuli)));
    angles = [angles; theseAngles(:)];
end
angles(:,2) = -1;
angles = reshape(angles',1,[])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(angles)

    stimulus = angles(j);
    if stimulus == -1
        ShowGreyScreen
    else
        showYourStimilus(stimulus);
    end
end

