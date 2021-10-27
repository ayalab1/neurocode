% [n, t] = ISIGrams(Res, Clu, SampleRate, BinSize, nBins)
%
% collect and plot a load of ISI histograms for several cells
%
% BinSize is in seconds
%
% n is a matrix n(TimeBin, CellNumber) giving the ISI count
% t is a time variable giving the time of the bins in seconds
%
% aza 2015, berenyi lab,
% From Harris et al 2002
% 
% Inputting Res = spikes.times (seconds)
%           Clu = spikes.UIDs

function [ISIs, n, t] = ISIGrams(Res, Clu, BinSize, nBins)

if (nargin<3) BinSize = 1e-3; end;
if (nargin<4) nBins = 100; end;

nCells = max(Clu);

n = zeros(nBins+1, nCells);

BinEdges = (0:nBins) * BinSize;
t = ((0:nBins-1) + 0.5)' * BinSize;
ISIs = cell(max(Clu), 1);

for Cell=1:nCells
    if length(Res{(Clu==Cell)})>=10 %minimum number of spikes
        ISIs{Cell} = diff(Res{(Clu==Cell)});
        n(:,Cell) = histc(ISIs{Cell}, BinEdges);
    else
        n(:,Cell) = zeros(length(BinEdges),1);
    end
end

% delete final bin because histc always produces a last bin with value 0
n(end,:) = [];
end