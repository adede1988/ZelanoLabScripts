function run_threshAnalysis_all(sessFilter, doKnit)
% RUN_THRESHANALYSIS_ALL  Master wrapper for the thresh-task macBP gamma analysis.
% =========================================================================
%   Runs the whole pipeline end-to-end, in dependency order, then knits the
%   R Markdown report. Ported from run_cueAnalysis_all; the one scientific change
%   is that the event-locked spectrograms and the time-resolved gamma are split by
%   ODOR CONDITION (behDat.type: low / med / air), all finalOnset-locked, instead
%   of by locking (trialStart / finalOnset).
%
%   USAGE
%     run_threshAnalysis_all()                 % full run over every thresh session
%     run_threshAnalysis_all({'<sessID>'})     % add/refresh ONE participant only
%     run_threshAnalysis_all([], false)        % full run, skip the R knit
%
%   sessFilter : [] = all on-disk thresh sessions; or a cellstr/string of sessIDs.
%   doKnit     : true (default) -> render threshTask_report.html via R at the end.
%
%   STAGES (each is an independent driver you can also run on its own)
%   (1) FOOOF / bestMac      -- run_thresh_tasks123 (all) OR run_thresh_fooof_one (one)
%   (2) noise + z-spectrograms (low/med/air) -- run_thresh_ztfr
%   (4) time-resolved gamma per condition    -- run_thresh_gamma_epochs
%   (3) group means + manifest (ALL sessions) -- run_thresh_task4_group
%   (5) report               -- threshTask_report.Rmd (knit in R)
%
%   WHERE RESULTS LIVE  (groupDir = THRESH_GROUPDIR or R:\...\Dupi_processing\groupStatFigs)
%     groupDir\threshTask_fooof_summary.csv   -- FOOOF + noise/QC (per channel)
%     groupDir\threshTask_gammaEpochs.csv     -- gamma peaks (session x cond x epoch)
%     groupDir\threshTask_data_manifest.csv   -- availability (derived)
%     groupDir\group_*_TFR_{low,med,air}.png, threshTask_group_means.mat
%     <id>\threshTask\  -- per-subject .png figures + <id>_thresh_bestMac_TFR.mat
%     bestMac is also written back INTO each <id>_PEA_threshold_preproc.mat
% =========================================================================

    if nargin < 1, sessFilter = []; end
    if nargin < 2 || isempty(doKnit), doKnit = true; end

    % self-bootstrap: make the repo (labPaths) + threshAnalysis reachable regardless
    % of the current folder, so thresh_init_paths()'s first labPaths() call resolves.
    here = fileparts(mfilename('fullpath'));      % ...\ZelanoLabScripts\threshAnalysis
    addpath(fileparts(here)); addpath(here);

    thresh_init_paths();              % paths once; the stages also call this (idempotent)
    L = labPaths();
    t0 = tic;

    % --- (1) FOOOF / bestMac ---------------------------------------------
    fprintf('\n########## (1) FOOOF / bestMac ##########\n');
    if isempty(sessFilter)
        run_thresh_tasks123([], true);
    else
        run_thresh_fooof_one(sessFilter, true);
    end

    % --- (2) noise + z-score spectrograms (low/med/air) ------------------
    fprintf('\n########## (2) noise + z-spectrograms ##########\n');
    run_thresh_ztfr(sessFilter);

    % --- (4) time-resolved gamma per condition (before the group step) ---
    fprintf('\n########## (4) time-resolved gamma ##########\n');
    run_thresh_gamma_epochs(sessFilter);

    % --- (3) group means + manifest (ALWAYS all sessions) ----------------
    fprintf('\n########## (3) group means + manifest ##########\n');
    run_thresh_task4_group([]);

    % --- (5) knit the report ---------------------------------------------
    rmd = fullfile(L.repo, 'threshAnalysis', 'threshTask_report.Rmd');
    if doKnit
        fprintf('\n########## (5) knit report ##########\n');
        knit_report(rmd);
    else
        fprintf('\n(5) skipped knit; render with:\n%s\n', knit_command(rmd));
    end

    fprintf('\nrun_threshAnalysis_all done in %.1f min.\n', toc(t0)/60);
end

% =========================================================================
function knit_report(rmd)
    cmd = knit_command(rmd);
    if isempty(cmd)
        fprintf(['  Rscript not found -- knit by hand in R:\n', ...
                 '    rmarkdown::render(''%s'', output_file=''threshTask_report.html'')\n'], ...
                 strrep(rmd,'\','/'));
        return;
    end
    fprintf('  %s\n', cmd);
    status = system(cmd);
    if status == 0
        fprintf('  report rendered: %s\n', fullfile(fileparts(rmd), 'threshTask_report.html'));
    else
        fprintf('  knit returned status %d -- run the command above manually if needed.\n', status);
    end
end

function cmd = knit_command(rmd)
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
            'rmarkdown::render(''', rmdF, ''', output_file=''threshTask_report.html'', quiet=TRUE); ', ...
            'cat(''RENDERED OK\n'')'];
    cmd = ['"' rscript '" -e "' expr '"'];
end
