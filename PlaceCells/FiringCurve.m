function [curve,stats] = FiringCurve(samples,spikes,varargin)

warning('FiringCurve is obsolete. Call FiringMap instead');
[curve,stats] = FiringMap(samples,spikes,varargin{:});


