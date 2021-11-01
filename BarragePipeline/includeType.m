function [keep] = includeType(ID, inc)
keep = zeros(length(ID),1); flg = 0;
for i = 1:length(ID)
    for j=1:length(inc)
        if ID(i) == inc(j)
            flg = flg+1;
        end
    end
    keep(i) = logical(flg);
    flg = 0;
end
end