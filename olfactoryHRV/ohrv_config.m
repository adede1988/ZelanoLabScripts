function P = ohrv_config()
%OHRV_CONFIG  Central paths for the olfactoryHRV analysis.
%
%   Every script in this folder calls this instead of hard-coding a path, so
%   the pipeline can be run from any machine or moved without edits.
%
%   P.work     intermediates: per-session slim extracts, beat files, CSVs
%   P.figs     generated figures
%   P.reports  generated HTML reports
%   P.dataRoot preprocessed session folders (the .mat finals)
%
%   To point at a different data root or scratch location, edit the two lines
%   below, or set the OHRV_WORK environment variable to override P.work.

here = fileparts(mfilename('fullpath'));

P.code     = here;
P.dataRoot = 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi';

w = getenv('OHRV_WORK');
if isempty(w), w = fullfile(here, 'work'); end
P.work    = w;
P.figs    = fullfile(here, 'figures');
P.reports = fullfile(here, 'reports');

for f = {'work','figs','reports'}
    if ~exist(P.(f{1}), 'dir'), mkdir(P.(f{1})); end
end
end
