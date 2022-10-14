function multiSaveFig(filename)
%multiSaveFig - saves figure as fig, png, and svg

set(gcf,'renderer','painters'); % otherwise svg can be saved as png if figure too big (uses OpenGL instead of painters)
[folder,name] = fileparts(filename);
filename = fullfile(folder,name); % gets rid of extension, if any
saveas(gcf,[filename '.fig']);
saveas(gcf,[filename '.png']);
saveas(gcf,[filename '.svg']);
% saveas(gcf,[filename '.pdf']); % These are almost always cut in a weird way; perform inkscape command to rpoduce pdf instead
try evalc(['!inkscape ' filename '.svg  --export-filename=' filename '.pdf']); end
end

