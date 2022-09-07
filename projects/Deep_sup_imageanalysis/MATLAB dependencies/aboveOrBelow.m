function [sign] = aboveOrBelow(pFoot,data_X,data_Y)
    %determine the sign of the depth. in CA1, sign==-1 indicates the cell
    %is deep. vice versa
    nCell = length(data_X);
    sign = ones(nCell,1);
    for i=1:nCell
        if(data_Y(i)<pFoot(i,2))
            sign(i) = -1;
        end
    end
end

