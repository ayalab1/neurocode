function [X1,X2,Y1,Y2] = expand_line(x1,x2,y1,y2)
        % Reflect p1 and p2 with respect to each other so that the section
        % between the new point P1 and P2 is larger (trippled.)
        X1 = x2-(x1-x2);
        X2 = x1-(x2-x1);
        Y1 = y2-(y1-y2);
        Y2 = y1-(y2-y1);
%         figure
%         hold on
%         plot([X1,X2],[Y1,Y2])
%         plot([x1,x2],[y1,y2])
end