function run_O15Analysis_all(sessFilter, doKnit)
% RUN_O15ANALYSIS_ALL  Master wrapper for the O15-task macBP gamma analysis.
% =========================================================================
%   Runs the whole pipeline end-to-end, in dependency order, then knits the
%   R Markdown report. The structure mirrors run_cueAnalysis_all, adapted to the
%   O15 task: events are SNIFFS split into three types (start / free / confirm),
%   and all spectrograms are z-scored against ONE shared baseline taken before
%   the first trial-start sniff (there is no cue-onset TTL to lock to in O15).
%
%   USAGE
%     run_O15Analysis_all()                 % full run over every FRESH O15 session
%     run_O15Analysis_all({'<sessID>'})     % add/refresh ONE participant only
%     run_O15Analysis_all([], false)        % full run, skip the R knit
%
%   sessFilter : [] = all fresh (modified >=2026-06-24) on-disk O15 sessions; or
%                a cellstr/string of sessIDs.
%   doKnit     : true (default) -> render O15Task_report.html via R at the end.
%
%   ONLY FRESH FINALS ARE USED. o15_session_table flags finals modified on/after
%   2026-06-24 (this O15 preprocessing run); older finals on disk are ignored.
%
%   STAGES (each is also an independent driver)
%   ---------------------------------------------------------------------
%   (1) FOOOF / bestMac          -- run_o15_fooof_all  OR  run_o15_fooof_one
%       FOOOFs every macBP channel (REUSES cue_fooof_macBP), picks bestMac (max
%       30-58 Hz periodic gamma peak; fallback: max flattened power), writes
%       bestMac INTO each <id>_O15preproc.mat, saves <id>_macBP_fooof_periodic.png,
%       and (re)writes O15Task_fooof_summary.csv. Sessions with no macBP are skipped.
%         Full (rewrites CSV):   run_o15_fooof_all([], true)
%         One new participant:   run_o15_fooof_one({'<sessID>'})
%
%   (2) Noise + baseline-z spectrograms -- run_o15_ztfr
%       Splits sniffs by sniffLabel, flags RELATIVE sharp-deflection noise
%       (cue_noise_trials, robust-z zd>K; K from o15_noise_K, tune via o15_calibrate_noise),
%       draws singleTrialRawMac_start.png, and computes bestMac baseline-z TF maps
%       for start/free/confirm (o15_ztfr_multi) vs the SINGLE pre-first-start
%       baseline, with respiration overlaid. Saves 3 PNGs + <id>_O15_bestMac_TFR.mat
%       and the noise/QC columns of O15Task_fooof_summary.csv.
%
%   (3) Time-resolved gamma      -- run_o15_gamma_epochs
%       bestMac, noise-clean, computed SEPARATELY for each sniff type (start/free/
%       confirm): Morlet power 5-58 Hz averaged into ten 250 ms windows
%       (-500..+2000 ms); FOOOF each; largest 30-58 Hz peak/window ->
%       O15Task_gammaEpochs.csv (col `locking`) + gammaTimeProgression_<type>.png.
%
%   (4) Group means + manifest   -- run_o15_task_group
%       Averages the per-subject z maps within each group (DupiS1/S2/S3/Control)
%       for the three lockings; saves group_<G>_O15TFR_<locking>.png +
%       O15Task_group_means.mat; rebuilds O15Task_data_manifest.csv. ALL sessions.
%
%   (5) Report                   -- O15Task_report.Rmd  (knit in R)
%
%   WHERE RESULTS LIVE  (R:\...\Lab_Common\Adam\Dupi_processing\)
%     groupStatFigs\O15Task_fooof_summary.csv   -- FOOOF + noise/QC (per channel)
%     groupStatFigs\O15Task_gammaEpochs.csv     -- gamma peaks (session x epoch)
%     groupStatFigs\O15Task_data_manifest.csv   -- availability (derived)
%     groupStatFigs\group_*_O15TFR_*.png, O15Task_group_means.mat
%     <id>\O15\  -- per-subject .png figures + <id>_O15_bestMac_TFR.mat
%     bestMac is also written back INTO each <id>_O15preproc.mat
% =========================================================================

    if nargin < 1, sessFilter = []; end
    if nargin < 2 || isempty(doKnit), doKnit = true; end

    here = fileparts(mfilename('fullpath'));      % ...\ZelanoLabScripts\O15Analysis
    addpath(fileparts(here)); addpath(here); addpath(fullfile(fileparts(here),'cueAnalysis'));

    o15_init_paths();
    L = labPaths();
    t0 = tic;

    fprintf('\n########## (1) FOOOF / bestMac ##########\n');
    if isempty(sessFilter)
        run_o15_fooof_all([], true);
    else
        run_o15_fooof_one(sessFilter, true);
    end

    fprintf('\n########## (2) noise + baseline-z spectrograms ##########\n');
    run_o15_ztfr(sessFilter);

    fprintf('\n########## (3) time-resolved gamma ##########\n');
    run_o15_gamma_epochs(sessFilter);

    fprintf('\n########## (4) group means + manifest ##########\n');
    run_o15_task_group([]);

    rmd = fullfile(L.repo, 'O15Analysis', 'O15Task_report.Rmd');
    if doKnit
        fprintf('\n########## (5) knit report ##########\n');
        knit_report(rmd);
    else
        fprintf('\n(5) skipped knit; render with:\n%s\n', knit_command(rmd));
    end

    fprintf('\nrun_O15Analysis_all done in %.1f min.\n', toc(t0)/60);
end

% =========================================================================
function knit_report(rmd)
    cmd = knit_command(rmd);
    if isempty(cmd)
        fprintf(['  Rscript not found -- knit by hand in R:\n', ...
                 '    rmarkdown::render(''%s'', output_file=''O15Task_report.html'')\n'], ...
                 strrep(rmd,'\','/'));
        return;
    end
    fprintf('  %s\n', cmd);
    status = system(cmd);
    if status == 0
        fprintf('  report rendered: %s\n', fullfile(fileparts(rmd), 'O15Task_report.html'));
    else
        fprintf('  knit returned status %d -- run the command above manually if needed.\n', status);
    end
end

function cmd = knit_command(rmd)
    cands = {'C:\Program Files\R\R-4.4.1\bin\Rscript.exe', ...
             'C:\Program Files\R\R-4.4.0\bin\Rscript.exe', 'Rscript'};
    rscript = '';
    for k = 1:numel(cands)
        if strcmp(cands{k},'Rscript') || isfile(cands{k}), rscript = cands{k}; break; end
    end
    if isempty(rscript), cmd = ''; return; end
    pandoc = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools';
    rmdF = strrep(rmd, '\', '/');
    expr = ['Sys.setenv(RSTUDIO_PANDOC=''', pandoc, '''); ', ...
            'rmarkdown::render(''', rmdF, ''', output_file=''O15Task_report.html'', quiet=TRUE); ', ...
            'cat(''RENDERED OK\n'')'];
    cmd = ['"' rscript '" -e "' expr '"'];
end
