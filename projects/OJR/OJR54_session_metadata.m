basepaths = {'Y:\OJRproject\OJR54\day6',...
    'Y:\OJRproject\OJR54\day8',...
    'Y:\OJRproject\OJR54\day10',...
    'Y:\OJRproject\OJR54\day11',...
    'Y:\OJRproject\OJR54\day12',...
    'Y:\OJRproject\OJR54\day13',...
    'Y:\OJRproject\OJR54\day14',...
    'Y:\OJRproject\OJR54\day15',...
    'Y:\OJRproject\OJR54\day16',...   
    'Y:\OJRproject\OJR54\day17',...
    'Y:\OJRproject\OJR54\day18',...
    };
 
 
% make basename.session
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.session.mat']));
%let the metatadata flow
    session.general.sessionType='Chronic';
    session.general.experimenters='Heath';
    session.general.location='Cornell';
    session.general.projects='OJR';
%brainRegions - CA1 fix extension on both  
    session.brainRegions.CA1.channels = [51,58,52,57,49,60,50,59,54,55,53,56,63,38,64,39,61,41,62,42,33,37,34,35,30,36,26,40,22,44,18,48,17,47,20,46,19,45,21,43,28,1,25,2,23,3,24,4,27,31,29,32,8,13,7,14,6,15,5,16,9,12,10,11];  
    session.brainRegions.CA1.electrodeGroups= [1,2,3,4,5];
    %PFC (PL)
    session.brainRegions.PL.channels=[65,78,83,96,79,82,66,80,81,77,84,68,76,85,69,75,86,70,74,87,95,73,88,94,72,89,93,71,90,92,91,67,97,110,115,128,111,114,98,112,113,99,109,116,100,108,117,101,107,118,102,106,119,127,105,120,126,104,121,125,103,122,124,123,129,142,147,160,143,146,130,144,145,131,141,148,132,140,149,133,139,150,134,138,151,159,137,152,158,136,153,157,135,154,156,155,161,174,179,192,175,178,162,176,177,163,173,180,164,172,181,165,171,182,166,170,183,191,169,184,190,168,185,189,167,186,188,187];
    session.brainRegions.PL.electrodeGroups=[6,7,8,9];
    %animal data    
    session.animal.name='OJR54';
    session.animal.sex='Male';
    session.animal.strain='C57Bl/6';
    session.animal.geneticLine='CamKII Cre:ChR2';
    %opticFibers
    session.animal.opticFiberImplants{1, 1}.opticFiber='105µm (125µm), NA 0.22 etched cone shaped (Thorlabs, FG105LCA)';
    session.animal.opticFiberImplants{1, 1}.supplier='Thorlabs';
    session.animal.opticFiberImplants{1, 1}.fiber_core_diameter=105;
    session.animal.opticFiberImplants{1, 1}.outer_diamter=125;
    session.animal.opticFiberImplants{1, 1}.numerical_aperature=0.2200;
    session.animal.opticFiberImplants{1, 1}.product_id=125;
    session.animal.opticFiberImplants{1, 1}.barinRegion='CA1';
    session.animal.opticFiberImplants{1, 1}.ap=2.3000;
    session.animal.opticFiberImplants{1, 1}.ml=1.6000;
    session.animal.opticFiberImplants{1, 1}.depth=1.1000;
    %surgery
    session.animal.surgeries{1, 1}.date='2022-06-10';
    session.animal.surgeries{1, 1}.start_time='10am';
    session.animal.surgeries{1, 1}.end_time='4pm';
    session.animal.surgeries{1, 1}.weight=35;
    session.animal.surgeries{1, 1}.type_of_surgery='Chronic';
    session.animal.surgeries{1, 1}.persons_involved='Heath';
    session.animal.surgeries{1, 1}.anesthesia='Heath';
    session.animal.surgeries{1, 1}.antibiotics='Heath';
    session.animal.surgeries{1, 1}.notes='Good, HPC somewhat inflammed';
    %probes -CA1
    session.animal.probeImplants{1, 1}.probe='A5x12-16-Buz-lin-5mm-100-200-160-177';
    session.animal.probeImplants{1, 1}.supplier='NeuroNexus';
    session.animal.probeImplants{1, 1}.nShanks=5;
    session.animal.probeImplants{1, 1}.nChannels=64;
    session.animal.probeImplants{1, 1}.shankSpacing=800;
    session.animal.probeImplants{1, 1}.layout='poly 2';
    session.animal.probeImplants{1, 1}.verticalSpacing='20';
    session.animal.probeImplants{1, 1}.descriptiveName='A5x12-16-Buz-lin-5mm-100-200-160-177 (64 ch, 5 shanks, poly 2 Non-uniform )';
    session.animal.probeImplants{1, 1}.brainRegion='CA1';
    session.animal.probeImplants{1, 1}.ap=2.3000;
    session.animal.probeImplants{1, 1}.ml=1.6000;
    session.animal.probeImplants{1, 1}.depth=1.1000;
    %probes - PFC (PL)
    session.animal.probeImplants{1, 2}.probe='PCB-128-5';
    session.animal.probeImplants{1, 2}.supplier='DiagnosticBiochips';
    session.animal.probeImplants{1, 2}.nShanks=4;
    session.animal.probeImplants{1, 2}.nChannels=128;
    session.animal.probeImplants{1, 2}.shankSpacing=150;
    session.animal.probeImplants{1, 2}.layout='poly 3';
    session.animal.probeImplants{1, 2}.verticalSpacing='20';
    session.animal.probeImplants{1, 2}.descriptiveName='PCB-128-5 (128 ch, 4 shanks, poly 3, 315um HD + staggered 600+200)';
    session.animal.probeImplants{1, 2}.brainRegion='PL';
    session.animal.probeImplants{1, 2}.ap=1.8000;
    session.animal.probeImplants{1, 2}.ml=0.3500;
    session.animal.probeImplants{1, 2}.depth=1.4000;
    %Tags (minimally used)
    session.analysisTags.probeslayout='poly2';
    session.analysisTags.probesVerticalSpacing=20;
   
   %epochs find certain epochs and add paradigms and envs
   for ii=1:length(session.epochs)
        if contains(session.epochs{1, ii}.name, "sleep")
            session.epochs{1, ii}.behavioralParadigm = 'Sleep';
            session.epochs{1, ii}.environment = 'Homecage';
        elseif contains(session.epochs{1, ii}.name, ["of", "openfield"])
            session.epochs{1, ii}.behavioralParadigm = 'Open Field';
            session.epochs{1, ii}.environment = 'Open Field';
        elseif contains(session.epochs{1, ii}.name, ["train" , "test"])
            session.epochs{1, ii}.behavioralParadigm = 'Spatial Object Recognition';
            session.epochs{1, ii}.environment = 'Open Field + Objects';
        elseif contains(session.epochs{1, ii}.name, "train")
            session.epochs{1, ii}.manipulation = 'SOR Training';
        elseif contains(session.epochs{1, ii}.name, "test")
            session.epochs{1, ii}.manipulation = 'SOR Test';
        else
            session.epochs{1, ii}.behavioralParadigm = '';
            session.epochs{1, ii}.environment = '';
        end
   end
    %spikeSorting
    session.spikeSorting{1, 1}.format='Phy';
    session.spikeSorting{1, 1}.method='KiloSort';
    session.spikeSorting{1, 1}.manuallyCurated=1;
   %session.spikeSorting{1, 1}.relativePath=''; implement to find kilosortfolder 
    %save updated session
    save(fullfile(basepath,[basename, '.session.mat']),'session');
   
end
%%
%Unique notes for individual animal
basepaths_ojr = {'Y:\OJRproject\OJR54\day6',...
    'Y:\OJRproject\OJR54\day8',...
    'Y:\OJRproject\OJR54\day10',...
    'Y:\OJRproject\OJR54\day11',...
    'Y:\OJRproject\OJR54\day12',...
    'Y:\OJRproject\OJR54\day13',...
    'Y:\OJRproject\OJR54\day14',...
    'Y:\OJRproject\OJR54\day16',...   
    'Y:\OJRproject\OJR54\day18',...
    };
 
man1="CL HPC SWR Extension";
man2="Random Delay HPC SWR Extension";
man3="No Stimulation";
stim1="";
stim2="CaMKII ChR2 LED";

manipulation=[man1,man3,man1,man1,man2,man1,man1,man1,man1];
stimulation=[stim2,stim1,stim2,stim2,stim2,stim2,stim2,stim2,stim2];
 
for hh=1:length(basepaths_ojr)
    basepath_ojr = basepaths_ojr{hh};
    basename_ojr = basenameFromBasepath(basepath_ojr);
    load(fullfile(basepath_ojr,[basename_ojr,'.session.mat']));
    for tt=1:length(session.epochs)
        if contains(session.epochs{1, tt}.name, "postsleep")
            session.epochs{1, tt}.manipulaton = manipulation(hh);
            session.epochs{1, tt}.stimuli = stimulation(hh);
        else
            session.epochs{1, tt}.behavioralParadigm = '';
            session.epochs{1, tt}.environment = '';
        end 
   end
end



