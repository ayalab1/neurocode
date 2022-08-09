function [] = darkBackground(hFigure,backColor,foreColor)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Summary: Quickly updates the color of a figure so that it can be easily
% displayed on a dark background, by changing background color and 
% converting axes, labels, etc to a light shade. However, if backColor and 
% foreColor are passed it can also be used to convert to other custom color
% schemes. By default (if none are specified in arguments), background is
% black and foreground is white.    
%
% Input Arguments: hFigure = valid handle for the figure to be converted
%                  backColor (optional) = RGB vector for desired background
%                                         color {black if none passed}
%                  foreColor (optional) = RGB vector for desired foreground
%                                         color {white if none passed}. To
%                                         set this, backColor must also be
%                                         passed to function.
% Output: none
%
% Origional Author: Elsa W. Birch, ewbirch@stanford.edu
%                   (feel free to modify to your own prefernces & share,
%                    but please give me credit for the base code, thanks!)
% Changes/Updates: 2011-01-19, inital coding
%
% More Info & Disclaimer:
%
% This may not be an entriely comprehensive tool, in fact it is most
% certainly not.  It is just my current quick and dirty solution to the
% fact that the energy barrier to my using pretty dark background (or any
% custom color) presentations is usually my plots.
%
% The default call (only input hFigure) gives somewhat similar results to 
% the command 'colordef black' and this function may in some respects be
% less comprehensive, but it works better for my own personal needs.
%
% Things it lacks color update for:
%   Lines/Patches/Markers - Basically I figure you desire, or at least
%                           have the capability, to set your own line
%                           colors. I know this function would not be
%                           useful for me if it overwrote all my line
%                           colors.  If you want an easy set of colors for
%                           a dark background use: 
%                           'set(gcf,'DefaultAxesColorOrder',colormap(summer(n)))'
%                           on your figure where n is the number of
%                           lineseries you have, or just remove argument to
%                           summer if you're using mesh/surf/contour type.
%                           A number of other predefined colormaps look
%                           good on a dark background.
%   Annotations - I assume if you are using arrows, etc you can just type
%                 in a color specification when you add them.
%   Plot Types - I have tested the function for some of the plot types I
%                most comonly use, however I certainly haven't gotten to
%                them all yet.  I will gradually, so check for updates
%
% Example 1: 
%
% % example data
% t = 0:0.1:10;
% y1 = sin(t);
% y2 = cos(t);
% % use plot
% f1 = figure;
% set(gcf,'DefaultAxesColorOrder',colormap(summer(2)))
% plot(t,y1,t,y2)
% L = legend('sine','cosine');
% title('testTitle')
% xlabel(' t must be time, right?')
% % Use dark grey & pink, just so its really ovbious what the function does
% darkBackground(f1,[0.2 0.2 0.2],[0.9 0.7 0.7])
% % Change the legend box, to demonstrate that you can modify any of the
% % objects changed by the darkBackground function after you run it.
% % (if you instead set 'box' 'off' for legend it is transparent)
% % set(L,'Edgecolor',pinkRGB)
% % Try commenting that line in/out just for fun, and if you want to see
% % somethign *really* fun, change the order of greyRGB and pinkRGB in the
% % call to darkBackground.
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Main function 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Start by checking that all arguments are of the correct type and set to
% default colors if none entered for backColor or foreColor
switch nargin
    case 0
        error('function requires arguments')
    case 1
        if strcmp(get(hFigure,'type'),'figure')
            % single input is figure handle, set colors to default
            backColor = [0 0 0]; % black
            foreColor = [1 1 1]; % white
        else
            error('argument is not valid figure handle')
        end
    case 2
        if strcmp(get(hFigure,'type'),'figure')
            if isValidRGB(backColor)
                % first input is figure handle, second is RGB, set
                % foreColor to default
                foreColor = [1 1 1]; % white
            else
                error('second argument must be valid RGB vector')
            end
        else
            error('first argument is not valid figure handle')
        end
    case 3
        if strcmp(get(hFigure,'type'),'figure')
            if isValidRGB(backColor)
                if isValidRGB(foreColor)
                    % all arguments valid, continue
                else
                   error('third argument must be valid RGB vector') 
                end
            else
                error('second argument must be valid RGB vector')
            end
        else
            error('first argument is not valid figure handle')
        end
    otherwise
        if strcmp(get(hFigure,'type'),'figure')
            if isValidRGB(backColor)
                if isValidRGB(foreColor)
                    % all arguments valid, continue
                    disp('warning: extra arguments ignored')
                else
                   error('third argument must be valid RGB vector') 
                end
            else
                error('second argument must be valid RGB vector')
            end
        else
            error('first argument is not valid figure handle')
        end
end % switch

% run function that recolors all relevent children in hierarchy of hFigure
recolorChildren(hFigure,backColor,foreColor)

end % main function darkBackground
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Sub functions 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% recolorChildren (Contains main functionality of darkBackground.m)
function [] = recolorChildren(hObject,backColor,foreColor)

% find all handles, including those that are hiddn of all children in the
% hierarchy under hObject
hChild = findall(hObject);

for iChild = 1:length(hChild);
    % get child type
    typeChild = get(hChild(iChild),'Type');
    % types of objects, and how to treat them
    switch typeChild
        case 'figure'
            % set figure background color
            set(hChild(iChild),'color',backColor)
        case 'axes'
            % Test axes tag since legends are also axes type, but require
            % differnt recoloring steps
            tagAxes = get(hChild(iChild),'Tag');
            if isempty(tagAxes)
                % Color: Color of the axes back planes. ({none} | ColorSpec)
                set(hChild(iChild),'Color','none')

                % XColor, YColor, ZColor: Color of axis lines, ticks, tick labels
                set(hChild(iChild),'XColor',foreColor)
                set(hChild(iChild),'YColor',foreColor)
                set(hChild(iChild),'ZColor',foreColor)

                % ? Visible: Axes/ticks visible ({on} | off)

                % Children:
                % image, light, line, patch, rectangle, surface, and text.
                % hidden(x-, y-, and z-axes and the title)
            elseif strcmp(tagAxes,'legend')
                % Color: Color of the axes back planes. ({none} | ColorSpec)
                % Make backColor in case of Legend, so that it is opaque if
                % it overlays plot elements
                set(hChild(iChild),'Color',backColor)
                % Although this seems like it should be taken care of by
                % text type, for a few legends it is not (no idea why)
                set(hChild(iChild),'TextColor',foreColor)
                % Set Edgecolor to foreColor
                set(hChild(iChild),'EdgeColor',foreColor)
            else
                % % % this is just a flag for while i am building this code
                % % disp(['other than empty or legend axes tag: ',tagAxes])
            end
        case 'surface'
            if ischar(get(hChild(iChild),'Facecolor'))
                % assume this means surface is a surf type plot, and leave
                % it alone
            else
                % assume in this case that surface is a mesh, and adjust
                % the facecolor to match background
                set(hChild(iChild),'FaceColor',backColor)
                % alpha adjustment here is a personal prefernce of mine
                set(hChild(iChild),'FaceAlpha',0.7)
            end
        case 'text'
                set(hChild(iChild),'Color',foreColor)
        otherwise
            % % % for function construction/debugging, display children not
            % % % recolored by function (but ignore user interaction types)
            % % if ~strcmp(typeChild(1:2),'ui')
            % % disp(['objects of Type = ',typeChild,...
            % %       ', not currently included in recolor'])
            % % end % if ~strcmp(typeChild(1:2),'ui')
    end % switch typeChild

end % for 

end % recolorChildren

% isValidRGB
function validRGB = isValidRGB(x)
% function tests if vector x is a valid 3 element RGB vector with values
% bounded by 0 and 1 (inclusive). Output is logical 
% 0 = not valid RGB
% 1 = valid RGB
sizeTest = all(size(x) == [1,3]);
valueTest = all((x >= 0) & (x <= 1));
if sizeTest && valueTest
    validRGB = 1;
else
    validRGB = 0;
end

end % isValidRGB