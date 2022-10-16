function animName = animalFromBasepath(basepath)
%% animalFromBasepath-Get animal name from basepath information
% Note that basepath should be the day folder of the animal that you're
% working with. The function will return the animal name as a string. 
% L Karaba, 10/21

remDay = find(basepath=='\');
cutBase = basepath(1:(remDay(end)-1)); %back out from day name
starAn = find(cutBase=='\');
aC = 1;
animName='';
for i = (starAn(end)+1):length(cutBase)
    animName(aC) = cutBase(i);
    aC = aC+1;
end
end