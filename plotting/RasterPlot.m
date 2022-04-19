function RasterPlot(spikes, height, varargin)

% "spikes" is a list of [timestamps id]
% provide desired height of spikes
% any other inputs will be passed on to "plot"

if nargin==1,
    height = 1;
end
    
times = spikes(:,1);
rows = spikes(:,2);

times = [times times nan(size(times))]';
rows =  [rows-0.45*height rows+0.45*height nan(size(rows))]';

plot(times(:),rows(:),varargin{:});