function [fig_width_pt, fig_height_pt] = set_size(width, fraction, subplots)
%     """Set figure dimensions to avoid scaling in LaTeX.
%     Parameters
%     ----------
%     width: float or string
%             Document width in points, or string of predined document type
%     fraction: float, optional
%             Fraction of the width which you wish the figure to occupy
%     subplots: array-like, optional
%             The number of rows and columns of subplots.
%     Returns
%     -------
%     fig_dim: tuple
%             Dimensions of figure in inches
%     """


if strcmp(width, 'thesis')
    width_pt = 426.79135;
elseif strcmp(width, 'beamer')
    width_pt = 307.28987;
elseif strcmp(width, 'paper'):
    width_pt = 595.276;
else
    width_pt = width;
end
%   Width of figure (in pts)
fig_width_pt = width_pt * fraction;

%   Golden ratio to set aesthetic figure height
%   https://disq.us/p/2940ij3
golden_ratio = (5^.5 - 1) / 2;

%   Figure height
fig_height_pt = fig_width_pt * golden_ratio * (subplots(1) / subplots(2));

end