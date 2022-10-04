% BEHAVIOR
%
% Files
%   align_tracking               - Heuristic to align behavior tracking trajectory with detected pulses
%   AngularVelocity              - AngularVelocity - Compute instantaneous angular velocity.
%   behaviorStructureDescription - AYA lab standard behavior structure
%   bz_scrubTracking             - USAGE  
%   check_for_dlc                - check_for_dlc - folder checker for dlc files
%   DefineZone                   - DefineZone - Define a restricted zone to which analyses can be circumscribed.
%   detectLED                    - Detect red and blue LEDs position in a video file and creates a 'led' file
%   find_and_move_videos         - matlab wrapper for find_and_move_videos.py
%   FindLapsNSMAadapted          - [laps, Vhs] = FindLaps_HorseShoe(V,newLapThreshold);
%   general_behavior_file        - converts multiple tracking data types to the standard described in cellexplorer
%   getSessionTracking           - Gets position trackign for each sub-session and concatenate all of them so they are
%   IsInZone                     - IsInZone - Test when the animal is in a given zone.
%   KalmanVel                    - root.KalmanVel(x, y, t, order);
%   LED2Tracking                 - Get LED tracking
%   linearTrackBehavior          - [behavior] = linearTrackBehavior(varargin)
%   LinearVelocity               - LinearVelocity - Compute instantaneous linear velocity.
%   load_dlc_csv                 - load_dlc_csv: general function to load dlc csv files
%   make_pre_task_post           - generate pre/task/post structure
%   manual_trackerjumps          - manual_trackerjumps: Allows you to manually cut out xy coordinates that
%   NSMAFindGoodLaps             - [startgoodlaps, stopgoodlaps, laps] = 
%   objPlaceScore                - [objScore] = objPlaceScore(sessionsSeq,basepath)
%   Plot_recording_States        - Plot recording states
%   process_and_sync_dlc         - Unpacks DLC CSV within subfolders
%   QuietPeriods                 - QuietPeriods - Find periods of immobility.
%   readOptitrackCSV             - INPUTS
%   RotateCoordinates            - angle should be in radians
%   socialPlusMazeBehavior       - [behavior] = socialPlusMazeBehavior(behavior,varargin)
%   Threshold                    - Threshold - Find periods above/below threshold.
%   trajectory_kalman_filter     - Kalman filter for obtaining an approximate Bayesian MAP estimate  to the rat's 
