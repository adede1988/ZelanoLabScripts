% run_chirp_job_lab.m  -- WMI-detached chirp pipeline launcher for lab E: session.
% Launch via:
%   ([wmiclass]'Win32_Process').Create('"matlab.exe" -batch "run(''E:/chirpAnalysis/run_chirp_job_lab.m'')"')
% Writes log to E:\chirp_job.log (sentinel E:\_CHIRP_DONE.txt on success).

diary('E:\chirp_job.log');
try
    fprintf('=== chirp job start %s ===\n', datestr(now));

    % 1. Code on path. The chirp code is the standalone E:\chirpAnalysis deploy, but
    %    setup_chirpAnalysis_paths needs labPaths() (lab repo) to resolve codePre=G: so
    %    it can add the external tools (Superlets/MPACT/FieldTrip). Put labPaths.m (lab
    %    repo) and labPaths_local.m (E:\ override -> codePre=G:, fieldtrip=E:) on the path
    %    FIRST, then E:\chirpAnalysis LAST so the freshest chirp_* win.
    addpath('G:\My Drive\GitHub\ZelanoLabScripts');   % labPaths.m + repo helpers
    addpath('E:\');                                   % labPaths_local.m (codePre=G:, fieldtrip=E:)
    addpath('E:\chirpAnalysis');                      % freshest chirp code (highest priority)
    setup_chirpAnalysis_paths(true);                  % now labPaths() resolves -> tools found
    assert(~isempty(which('mp_adapt_chirplets')), 'MPACT not on path -- chirplet test would be empty');

    % 2. E:-local data + output dirs. FIGROOT mirrors the home subject-fig layout
    %    (<root>\Adam\Dupi_processing\<id>\) so figs copy back to R: cleanly; each
    %    session's single-trial TF plots land in <id>\singleTrialTF\ (C.figSub).
    setenv('CHIRP_DATAROOT', 'E:\Lab_Common');
    setenv('CHIRP_OUTDIR',   'E:\chirpOut');
    setenv('CHIRP_FIGROOT',  'E:\Lab_Common\Adam\Dupi_processing');

    % 2b. Clear stale per-trial CSVs (they are APPENDED across the batch -> a re-run would
    %     otherwise double-count). subject_summary.csv is overwritten, but clear it too.
    od = getenv('CHIRP_OUTDIR');
    if ~isempty(od) && isfolder(od)
        stale = [dir(fullfile(od,'*_trials.csv')); dir(fullfile(od,'subject_summary.csv')); ...
                 dir(fullfile(od,'verdict_vectors.csv'))];
        for k = 1:numel(stale)
            try, delete(fullfile(od, stale(k).name)); catch, end
        end
        fprintf('cleared %d stale CSV(s) in %s\n', numel(stale), od);
    end

    % 3. Run pipeline on all fresh cue sessions (freshOnly=true, saveFigs=true)
    maxNumCompThreads(6);   % leave 2 threads free
    opts.saveFigs = true;
    opts.nFigTrials = 10;   % per-session: up to 10 single-trial spectrograms
    run_chirp_analysis_pipeline([], opts);

    % 4. Group adjudication (phase 2): per-session verdict vectors + group tests.
    try
        chirp_group_adjudicate(od, chirp_config());
    catch GE
        fprintf('group adjudication failed: %s\n', GE.message);
    end

    fprintf('=== chirp job DONE %s ===\n', datestr(now));
    fid = fopen('E:\_CHIRP_DONE.txt','w'); fprintf(fid,'done\n'); fclose(fid);
catch ME
    fprintf('ERROR: %s\n%s\n', ME.message, getReport(ME,'extended'));
    fid = fopen('E:\_CHIRP_ERROR.txt','w'); fprintf(fid,'%s\n',ME.message); fclose(fid);
end
diary off
