function nNew = run_o15_integrate_new(doKnit)
% RUN_O15_INTEGRATE_NEW  Detect O15 finals that are fresh on disk but NOT yet in
%   the analysis, and integrate each via the single-session path. Used by the
%   every-10-minute monitoring loop while O15 preprocessing is still running.
%
%   nNew = run_o15_integrate_new(doKnit)
%     doKnit : re-render O15Task_report.html after integrating (default true).
%
%   "New" = a fresh (modified >=2026-06-24) on-disk session whose sessID is not
%   already present in O15Task_data_manifest.csv. Each new session is run through
%   run_O15Analysis_all({id}) -> FOOOF(append) + baseline-z spectrograms + gamma
%   + group means/manifest(all) + (optional) knit. Returns the count integrated
%   (0 = nothing new this tick).

    if nargin < 1 || isempty(doKnit), doKnit = true; end
    % LIGHTWEIGHT detection path: only addpath the repo (NOT eeglab/fieldtrip) so
    % a "nothing new" tick is cheap. EEGLAB/FieldTrip are initialised by
    % run_O15Analysis_all only when there is actually a new session to process.
    here = fileparts(mfilename('fullpath')); repo = fileparts(here);
    addpath(repo); addpath(fullfile(repo,'cueAnalysis')); addpath(here);
    groupDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs';
    manPath  = fullfile(groupDir, 'O15Task_data_manifest.csv');

    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);
    done = strings(0,1);
    if isfile(manPath)
        try, M = readtable(manPath); done = string(M.sessID); catch, end
    end
    newIDs = setdiff(T.sessID, done);

    if isempty(newIDs)
        nNew = 0;
        fprintf('[%s] no new fresh O15 finals (%d fresh on disk, all integrated).\n', ...
            datestr(now,'yyyy-mm-dd HH:MM'), height(T));
        return;
    end
    nNew = numel(newIDs);
    fprintf('[%s] integrating %d new O15 session(s): %s\n', ...
        datestr(now,'yyyy-mm-dd HH:MM'), nNew, strjoin(cellstr(newIDs),', '));

    run_O15Analysis_all(cellstr(newIDs), doKnit);
    fprintf('[%s] integration of %d session(s) complete.\n', datestr(now,'yyyy-mm-dd HH:MM'), nNew);
end
