function PlotCircle(x,y,r,varargin)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
% https://fr.mathworks.com/matlabcentral/answers/3058-plotting-circles
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,varargin{:});
end