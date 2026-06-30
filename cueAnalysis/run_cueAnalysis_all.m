function run_cueAnalysis_all(sessFilter, doKnit)
% RUN_CUEANALYSIS_ALL  Master wrapper for the cue-task macBP gamma analysis.
% =========================================================================
%   Runs the whole pipeline end-to-end, in dependency order, then knits the
%   R Markdown report. There is otherwise NO single entry point -- each stage
%   below is an independent driver you can also run on its own.
%
%   USAGE
%     run_cueAnalysis_all()                 % full run over every cue session
%     run_cueAnalysis_all({'<sessID>'})     % add/refresh ONE participant only
%     run_cueAnalysis_all([], false)        % full run, skip the R knit
%
%   sessFilter : [] = all on-disk cue sessions; or a cellstr/string of sessIDs.
%   doKnit     : true (default) -> render cueTask_report.html via R at the end.
%
%   WHAT EACH STAGE DOES, AND HOW TO RUN IT INDEPENDENTLY
%   ---------------------------------------------------------------------
%   (1) FOOOF / bestMac          -- run_cue_tasks123  OR  run_cue_fooof_one
%       FOOOFs every macBP channel (FieldTrip ft_freqanalysis -> process_fooof),
%       picks the channel with the largest 30-58 Hz periodic gamma peak as
%       `bestMac` (fallback: largest flattened power if no channel has a peak),
%       writes `bestMac` INTO each <id>_cueTaskPreproc.mat, saves the per-subject
%       <id>_macBP_fooof_periodic.png, and (re)writes cueTask_fooof_summary.csv.
%         Standalone (all sessions, REWRITES the whole CSV):
%             run_cue_tasks123([], true)
%         Standalone (ONE new participant, APPENDS without touching others):
%             run_cue_fooof_one({'260514_OBE_NWU_BW_1'})
%       --> MUST run before stages 2 and 4 (they read the saved `bestMac`).
%
%   (2) Noise + z-score spectrograms -- run_cue_ztfr
%       Flags noisy trials (sharp-deflection rule: >80 uV max-min in any 10 ms
%       window of the trialStart-locked bestMac epoch -> cue_noise_trials), draws
%       the red-flagged singleTrialRawMac.png, then computes per-trial bootstrap-z
%       (myChanZscore) time-frequency maps (EEGLAB newtimef) of bestMac locked to
%       trialStart and to finalOnset (both vs the SAME trialStart baseline), with
%       the mean respiration overlaid. Saves <id>_<bestMac>_TFR_{trialStart,
%       finalOnset}.png + the numeric <id>_cue_bestMac_TFR.mat, and writes the
%       noise/QC columns onto the bestMac row of cueTask_fooof_summary.csv.
%         Standalone:  run_cue_ztfr([])   % or run_cue_ztfr({'<sessID>'})
%       --> needs stage 1 first; MUST run before stage 3 (it averages its .mat).
%
%   (3) Group-mean spectrograms  -- run_cue_task4_group
%       Loads every per-subject <id>_cue_bestMac_TFR.mat, averages the z-maps and
%       respiration within each group (DupiS1/S2/S3/Control) on a shared +-5 z
%       scale, saves the 8 group_<Group>_TFR_<locking>.png + cueTask_group_means
%       .mat, and rebuilds the data-availability manifest (cue_make_manifest).
%       ALWAYS aggregates ALL sessions (ignore sessFilter here).
%         Standalone:  run_cue_task4_group([])
%       --> needs stage 2 (the per-subject TFR .mat files) to exist.
%
%   (4) Time-resolved gamma      -- run_cue_gamma_epochs
%       For bestMac on the noise-clean finalOnset-locked trials: Morlet power
%       5-58 Hz (newtimef), averaged into ten 250 ms windows (-500..+2000 ms);
%       FOOOFs each window; records the largest 30-58 Hz peak per window to
%       cueTask_gammaEpochs.csv and saves the ochre->purple gammaTimeProgression
%       .png. Independent of stages 2-3 (only needs `bestMac` from stage 1).
%         Standalone:  run_cue_gamma_epochs([])   % or {'<sessID>'}
%
%   (5) Report                   -- cueTask_report.Rmd  (knit in R)
%       Reads the three CSVs + the saved figures and renders cueTask_report.html
%       (methods, data table, FOOOF table, noise table + example plots, single-
%       subject and group spectrograms, and the time-resolved gamma section with
%       a per-group linear mixed-effects test). Knit by this wrapper (or by hand,
%       see knit_report below).
%
%   UTILITY SCRIPTS (not part of the linear run, but handy)
%     run_cue_singletrial  -- re-plot ONLY singleTrialRawMac.png (no spectrograms)
%     run_cue_fooof_one    -- FOOOF + append one session (used by stage 1 above)
%     cue_make_manifest    -- regenerate cueTask_data_manifest.csv from the CSV
%
%   SHARED HELPERS (called by the drivers; rarely run directly)
%     cue_init_paths     -- EEGLAB + FieldTrip(+brainstorm FOOOF) onto the path
%     cue_session_table  -- session list/group/paths from dataTracking.xlsx
%     cue_fooof_macBP    -- FOOOF a session's macBP channels -> bestMac
%     cue_noise_trials   -- per-trial sharp-deflection noise flags
%     cue_ztfr_pair      -- the bootstrap-z TFR engine (both lockings)
%     cue_plot_{fooof,ztfr,singletrial}  -- the figure plotters
%     (repo-level: labPaths, applyParams, myChanZscore)
%
%   WHERE RESULTS LIVE  (R:\...\Lab_Common\Adam\Dupi_processing\)
%     groupStatFigs\cueTask_fooof_summary.csv   -- FOOOF + noise/QC (per channel)
%     groupStatFigs\cueTask_gammaEpochs.csv     -- gamma peaks (session x epoch)
%     groupStatFigs\cueTask_data_manifest.csv   -- availability (derived)
%     groupStatFigs\group_*_TFR_*.png, *_group_means.mat
%     <id>\cueTask\  -- per-subject .png figures + <id>_cue_bestMac_TFR.mat
%     bestMac is also written back INTO each <id>_cueTaskPreproc.mat
% =========================================================================

    if nargin < 1, sessFilter = []; end
    if nargin < 2 || isempty(doKnit), doKnit = true; end

    % self-bootstrap: make the repo (labPaths) + cueAnalysis reachable regardless
    % of the current folder, so cue_init_paths()'s first labPaths() call resolves.
    here = fileparts(mfilename('fullpath'));      % ...\ZelanoLabScripts\cueAnalysis
    addpath(fileparts(here)); addpath(here);

    cue_init_paths();                 % paths once; the stages also call this (idempotent)
    L = labPaths();
    t0 = tic;

    % --- (1) FOOOF / bestMac ---------------------------------------------
    % Full run -> run_cue_tasks123 (rewrites the whole CSV). Single/subset run
    % -> run_cue_fooof_one (appends, so the other sessions' rows are preserved).
    fprintf('\n########## (1) FOOOF / bestMac ##########\n');
    if isempty(sessFilter)
        run_cue_tasks123([], true);
    else
        run_cue_fooof_one(sessFilter, true);
    end

    % --- (2) noise + z-score spectrograms --------------------------------
    fprintf('\n########## (2) noise + z-spectrograms ##########\n');
    run_cue_ztfr(sessFilter);

    % --- (4) time-resolved gamma (independent of 2/3; do before the group step) --
    fprintf('\n########## (4) time-resolved gamma ##########\n');
    run_cue_gamma_epochs(sessFilter);

    % --- (3) group means + manifest (ALWAYS all sessions) ----------------
    fprintf('\n########## (3) group means + manifest ##########\n');
    run_cue_task4_group([]);

    % --- (5) knit the report ---------------------------------------------
    rmd = fullfile(L.repo, 'cueAnalysis', 'cueTask_report.Rmd');
    if doKnit
        fprintf('\n########## (5) knit report ##########\n');
        knit_report(rmd);
    else
        fprintf('\n(5) skipped knit; render with:\n%s\n', knit_command(rmd));
    end

    fprintf('\nrun_cueAnalysis_all done in %.1f min.\n', toc(t0)/60);
end

% =========================================================================
function knit_report(rmd)
% Render the Rmd to HTML via Rscript. Tries common Rscript locations; if none
% is found (or it errors), prints the exact command to run by hand.
    cmd = knit_command(rmd);
    if isempty(cmd)
        fprintf(['  Rscript not found -- knit by hand in R:\n', ...
                 '    rmarkdown::render(''%s'', output_file=''cueTask_report.html'')\n'], ...
                 strrep(rmd,'\','/'));
        return;
    end
    fprintf('  %s\n', cmd);
    status = system(cmd);
    if status == 0
        fprintf('  report rendered: %s\n', fullfile(fileparts(rmd), 'cueTask_report.html'));
    else
        fprintf('  knit returned status %d -- run the command above manually if needed.\n', status);
    end
end

function cmd = knit_command(rmd)
% Build the Rscript render command (or '' if Rscript can't be located).
    cands = {'C:\Program Files\R\R-4.5.2\bin\Rscript.exe', ...
             'C:\Program Files\R\R-4.4.1\bin\Rscript.exe', ...
             'C:\Program Files\R\R-4.4.0\bin\Rscript.exe', 'Rscript'};
    rscript = '';
    for k = 1:numel(cands)
        if strcmp(cands{k},'Rscript') || isfile(cands{k}), rscript = cands{k}; break; end
    end
    if isempty(rscript), cmd = ''; return; end
    pandoc = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools';
    rmdF = strrep(rmd, '\', '/');
    expr = ['Sys.setenv(RSTUDIO_PANDOC=''', pandoc, '''); ', ...
            'rmarkdown::render(''', rmdF, ''', output_file=''cueTask_report.html'', quiet=TRUE); ', ...
            'cat(''RENDERED OK\n'')'];
    cmd = ['"' rscript '" -e "' expr '"'];
end
