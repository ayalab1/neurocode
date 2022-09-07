function h = yline(vals)
% create horizontal lines at vals

    hold on
    h = {};
    lims = get(gca,'xlim');
    for kp = 1:length(vals)
        h{kp} = plot(lims, [vals(kp) vals(kp)], '--k');
    end
    set(gca,'xlim', lims)

end