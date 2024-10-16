function [fig_width_px, fig_height_px] = set_size(width, fraction, subplots)
%     """Set figure dimensions to avoid scaling in LaTeX.
%     Parameters
%     ----------
%     width: float or string
%             Document width in pixels, or string of predined document type
%     fraction: float, optional
%             Fraction of the width which you wish the figure to occupy
%     subplots: array-like, optional
%             The number of rows and columns of subplots.
%     Returns
%     -------
%     fig_dim: tuple
%             Dimensions of figure in px
%     """


if strcmp(width, 'thesis')
    width_px = 569.055;
elseif strcmp(width, 'beamer')
    width_px = 409.720;
elseif strcmp(width, 'paper')
    width_px = 793.701;
else
    width_px = width;
end
%   Width of figure (in pts)
fig_width_px = width_px * fraction;

%   Golden ratio to set aesthetic figure height
%   https://disq.us/p/2940ij3
golden_ratio = (5^.5 - 1) / 2;

%   Figure height
fig_height_px = fig_width_px * golden_ratio * (subplots(1) / subplots(2));

end