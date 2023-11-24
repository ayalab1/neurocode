function [sequence,shuffled] = ReplayScore(matrix,varargin)

% ReplayScore
% This is a structure wrapper which calls "FindReplayScore"
%
% USAGE
%
%    [r,p,a,b,rShuffled,c,cShuffled,jump,jumpShuffled,maxJump, maxJumpShuffled] = ReplayScore(matrix,threshold,<options>);
%
%  INPUT
%
%    matrix     probability matrix for a specific event ("estimations" output of ReconstructPosition)       
%    threshold  distance away from the fitted line (in bins) to count towards the score (see Davidson et al. 2009)
%    
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   considered distance from the line (default = 15 bins)
%     'nShuffles'   default = 500
%     'circular'    for circular-linear data (default = 'off')
%     'shuffle'     either 'column' (spatial) or 'temporal' (default = 'column');
%    =========================================================================
%
%   OUTPUT
%
%     r               replay score of matrix
%     p               p-value of replay score (based on a column shuffle)
%     start           start position bin of the fitted line
%     stop            stop position bin of the fitted line
%     slope           slope of the fitted line (in bin units)
%     rShuffled       replay scores of shuffled matrices
%     stShuffled      st of shuffled matrices
%     spShuffled      sp of shuffled matrices
%     wc              weigthed correlation of matrix
%     wcShuffled       weigthed correlations of shuffled matrices
%     jump            jump value of matrix  
%     jumpShuffled    jump values of shuffled matrices
%     maxJump         max jump value of matrix
%     maxJumpShuffled max jump values of shuffled matrices
%
%  SEE ALSO
%
% Copyright (C) 2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%------------------------------------------------------------------------

[r,p,st,sp,rShuffled,aShuffled,bShuffled,c,cShuffled,jump,jumpShuffled,maxJump,maxJumpShuffled,quadrantScore,qShuffled] = FindReplayScore(matrix,varargin{:});

sequence.quadrantScore = quadrantScore;
sequence.score = r;
sequence.zscore = (r - nanmean(rShuffled))./nanstd(rShuffled);
sequence.pValue = p; sequence.pValue(isnan(r)) = nan;
sequence.lineStart = st;
sequence.lineStop = sp;
sequence.slope = (sp-st)./size(matrix,2);
sequence.weightedCorrelation = c;
sequence.jump.mean = jump;
sequence.jump.max = maxJump;
sequence.zWeighted = (abs(c) - mean(abs(cShuffled),2))./std(abs(cShuffled),[],2);
sequence.pWeighted = sum(abs(cShuffled)>=abs(c),2)./sum(~isnan(cShuffled),2);
sequence.pWeighted(isnan(c) | nanstd(abs(cShuffled))==0) = nan;
sequence.zSignedWeighted = (c - mean(cShuffled,2))./std(cShuffled,[],2);
sequence.pSignedWeighted = sum(cShuffled>=abs(c),2)./sum(~isnan(cShuffled),2);
sequence.pSignedWeighted(isnan(c) | nanstd(cShuffled)==0) = nan;
sequence.jump.zMax = (maxJump - mean(maxJumpShuffled,2))./std(maxJumpShuffled,[],2);
sequence.jump.pMax = sum(maxJumpShuffled<sequence.jump.max,2)/size(maxJumpShuffled,2);
sequence.jump.zMean = (jump - mean(jumpShuffled,2))./std(jumpShuffled,[],2);
sequence.jump.pMean = sum(jumpShuffled < jump,2)/size(jumpShuffled,2);
sequence.pQuadrant =  sum(qShuffled>=quadrantScore,2)./sum(~isnan(qShuffled),2);
shuffled.score = rShuffled;
shuffled.lineStart = aShuffled;
shuffled.lineStop = bShuffled;
shuffled.weightedCorrelation = cShuffled;
shuffled.zWeighted = cShuffled;
shuffled.jump.mean = jumpShuffled;
shuffled.jump.max = maxJumpShuffled;







