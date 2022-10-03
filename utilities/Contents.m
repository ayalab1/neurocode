% UTILITIES
%
% Files
%   Accumulate                 - Accumulate - Accumulate repeated observations.
%   animalFromBasepath         - Get animal name from basepath information
%   basenameFromBasepath       - gets basename from a basepath input to a function.  Uses name of the last
%   Bin                        - Bin - Assign each input value a bin number between 1 and N.
%   Bins                       - Bins - build temporal bins with specifying start and end points, window and step.  
%   bz_BimodalThresh           - [thresh,cross,bihist,diptest,overthresh] = bz_BimodalThresh(bimodaldata) 
%   bz_Counter                 - bz_Counter(current,total,name) counts up loop indices, given current index
%   bz_Diff                    - Diff - Differentiate.
%   bz_hartigansdipsigniftest  - function		[dip,p_value,xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
%   bz_hartigansdiptest        - USAGE
%   bz_IDXtoINT                - bz_IDXtoINT(IDX) Converts state indices to state on/offsets
%   bz_INTtoIDX                - [IDX] = bz_INTtoIDX(INT) Converts state on/offsets to vector of indices
%   bz_isSessionInfo           - USAGE
%   bz_NormToRange             - [ normdata ] = bz_NormToRange(data,range,databounds) normalizes some data to
%   bz_RandomWindowInIntervals - [ window ] = bz_RandomWindowInIntervals( interval,winsize,nwins ) returns a random
%   bz_SpktToSpkmat            - spikemat = bz_SpktToSpkmat(spiketimes,<options>) takes a 
%   CatCon                     - Concatenation with Condition function
%   catOverlap                 - catOverlap - concatenate overlapping timebins
%   checkFile                  - Checks that file exists and returns proper basepath and filename
%   CircularShift              - CircularShift - Shift matrix rows or columns circularly.
%   Clip                       - Clip - Clip values.
%   collapseStruct             - structout = CollapseStruct( structin,dim,combine,NEST ) Combines elements in a
%   Concentration              - Concentration - Estimate the concentration parameter for circular data.
%   ConsolidateIntervalsFast   - ConsolidateIntervalsFast - a faster version of ConsolidateIntervals.
%   CumSum                     - CumSum - Cumulative sum of elements. Partial sums can also be computed.
%   Diff                       - Diff - Differentiate.
%   Dist                       - Outputs the normalised distributions of the data provided at a given resolution.
%   eventIntervals             - EVENTINT Separate given  event struct by given time intervals
%   fastrms                    - Instantaneous root-mean-square (RMS) power via convolution.
%   FConv                      - FConv(kernel,signal) convolves a signal with a kernal via the fourier
%   FindClosest                - FindClosest - field the indices of a reference vector closest to a query vector
%   FindFieldHelper            - [field] =  FindFieldHelper(map,x,y,threshold,circX,circY)
%   findIntervals              - 
%   FindIntsNextToInts         - [ints1,ints2,ints1idx,ints2idx] = FindIntsNextToInts(ints1,ints2,tol) 
%   FindLocalMaxima            - 
%   FindLocalMinima            - 
%   findmax                    - 
%   findmin                    - 
%   gausskernel                - create a 2D gaussian.
%   getCurState                - comp = ['WAKEstate', 'NREMstate', 'REMstate', 'THETA', 'nonTHETA', 'THETAtask'];
%   GLMgain                    - GLMgain - Train a GLM on your "source" data and see how well it predicts 
%   Group                      - Provide as many vectors as you want. They will be grouped in a matrix
%   inRange                    - 
%   Insert                     - Insert - Insert lines in a matrix.
%   Interpolate                - Interpolate - Interpolate samples (positions, spikes, LFP, etc.) at given timestamps.
%   intervalsBetween           - intervalsBetween - Pull times between event periods
%   ISIGrams                   - [n, t] = ISIGrams(Res, Clu, SampleRate, BinSize, nBins)
%   JointPETH                  - JointPETH - produce a joint histogram for the co-occurrence of two sets of signals around events. 
%   linspaceVector             - linspaceVector - produce linearly spaced vector with multiple [start stop] points. 
%   LoadXml                    - function [xml, rxml] = LoadXml(FileBase)
%   Match                      - Match - Replace values in one list with closest values in a second list.
%   matchfields                - [ newstruct1, newstruct2 ] = bz_Matchfields( struct1,struct2,mode )
%   mean2str                   - mean2str - Convert mean or median and SEM (or confidence interval) to string.
%   mirror2NaN                 - INPUTS
%   mmax                       - mmax - Return overall min and max values for any number of arrays.
%   nancorr                    - nancorr - call <a href="matlab:help corr">corr</a> ignoring nans
%   nansem                     - nansem - Compute standard error of the mean (SEM), ignoring NaNs.
%   nansmooth                  - Smooth - Smooth using a Gaussian kernel.
%   nanzscore                  - zscore that ignores nans
%   NormToInt                  - NormToInt(data,normtype, int,sf) normalizes the data to a subset of time 
%   out2                       - returns second output argument of fun
%   parfor_wait                - This class creates a waitbar or message when using for or parfor.
%   PlotColorMap               - PlotColorMap - Plot a color map.
%   Portion                    - 
%   RandomIntervals            - RandomIntervals - generates n random intervals from 0 to 1, covering a total of 'sumDuration'
%   RelativeTimes              - Relative times - compute the relative times of samples within intervals
%   RemoveOutliers             - 
%   RemoveRepetitions          - RemoveRepetitions - Removes repetitions from a vector or a matrix.
%   Renumber                   - This function takes a vector and renumbers all the entries so there are
%   rep_zero                   - 
%   Restrict                   - Restrict - Keep only samples that fall in a given list of time intervals.
%   Rsquared                   - Rsquared - return the r^2 value of a fit between a linear model and data
%   sec2min                    - 
%   sem                        - sem - Compute standard error of the mean (SEM).
%   semedian                   - semedian - Compute standard error of the median.
%   Shrink                     - Shrink - reduce matrix resolution by shrinking it (akin to pixelating)
%   Shuffle                    - Shuffle - shuffle the elements of a matrix
%   Smooth                     - Smooth - Smooth using a Gaussian kernel.
%   SpikeSuprise               - SpikeSuprise - estimate the surprise of a spike count assuming poisson spiking
%   Sync                       - Sync - Make sample timestamps relative to synchronizing events.
%   SyncHist                   - SyncHist - Compute a histogram on event-synchronized samples (e.g. a PSTH).
%   Unfind                     - Unfind - the inverse operation of 'find': x = find(y); y = Unfind(x);
%   Unshift                    - The opposite of the option 'shift' in Restrict
%   v2struct                   - v2struct
%   xmltools                   - - tools for managing xml data sets
%   ZeroCrossings              - ZeroCrossings - Test zero crossings in a given time series.
%   ZeroToOne                  - ZeroToOne - Normalize values in [0,1].
