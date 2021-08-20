%% Basic pipeline for place field analysis 
    % These are the basic fucntions for performing place field analysis (in
    % 1D for now). Is still work in progress and anyone is welcome to contribute. 
    
    % TODO:
    % - Decide a common trial structure and way to linearize positions.  
    % - Implent trial based firing map and place field detection
    % - 2D place fields
    % - phase precession 
    % - theta compression and sequences
    % - position decoding 

   % 0- Linearize maze positions and get trial structure
   
   % 1-Calculate firing maps 
    [firingMaps] = firingMapAvg(posTrials,spikes);
    save([basename '.firingMapsAvg.cellinfo.mat'],'firingMaps');
    
   % 2- Detect (average) place fields
    [placeFieldStats] = findPlaceFieldsAvg1D('firingMaps',firingMaps,'minPeak',1,'sepEdge',0.04);
    save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');
    
   % 3- Calculate template of place fields in the maze (for decoding, etc.)
    [placeFieldTemplate] = findPlaceFieldsTemplate('placeFieldStats',placeFieldStats,'firingMaps',firingMaps);
    save([basename '.placeFieldTemplate.mat'],'placeFieldTemplate');
   

