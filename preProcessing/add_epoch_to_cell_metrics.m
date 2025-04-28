function add_epoch_to_cell_metrics(basepath)
% add_epoch_to_cell_metrics: adds epochs to cell_metrics for viewing in gui
% for case when cell_metrics was created without access to session file

basename = basenameFromBasepath(basepath);

load(fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']), 'cell_metrics')
load(fullfile(basepath, [basename, '.session.mat']), 'session')


cell_metrics.general.epochs = session.epochs;
save(fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']), 'cell_metrics')

end
