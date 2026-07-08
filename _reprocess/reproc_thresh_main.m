function reproc_thresh_main(idxList)
% REPROC_THRESH_MAIN  Portable main step (intermediate -> final) for threshTask,
% mirroring threshPreProc_main.m (incl. the 45-trial sniff-TTL rebuild). Run
% AFTER reproc_backup('threshTask') + reproc_thresh_make. Curated-only, resumes
% via done-marker, per-session failure isolation (incl. ZLP:blinkAmbiguous).
%
%   reproc_thresh_main()         % all curated thresh sessions
%   reproc_thresh_main(idxList)  % only these indices

    L = labPaths();
    addpath(genpath(L.repo));
    addpath(genpath(L.eeglab));
    figPath = L.figPath;
    EEGLOC  = readtable(L.eegLocCsv);
    set(0, 'defaultfigurevisible', 'off');

    logf     = fopen(fullfile(reproc_root(), 'thresh_run.log'), 'a');
    donePath = fullfile(reproc_root(), 'thresh_done.txt');
    doneIds  = reproc_readDone(donePath);

    cfg = applyParams('threshTask', 'main');
    if nargin < 1 || isempty(idxList), idxList = find(~strcmpi(cfg.paramSource, 'guess')); end
    logmsg(logf, sprintf('==== reproc_thresh_main start %s ====', datestr(now)));

    for s = idxList(:)'
        id = cfg.sessionIDs{s};
        if strcmpi(cfg.paramSource{s}, 'guess'), logmsg(logf, sprintf('GUESS_SKIP_thresh %s', id)); continue; end
        if any(strcmpi(id, doneIds)),            logmsg(logf, sprintf('ALREADY_DONE_thresh %s', id)); continue; end
        try
            S = struct('id', id, 'root', cfg.root{s}, 'fig', fullfile(figPath, id));
            preDir = fullfile(S.root, S.id, 'preProc');
            matf   = fullfile(preDir, [S.id '_PEA_threshold_preproc.mat']);
            if ~exist(matf, 'file'), logmsg(logf, sprintf('NOINT_thresh %s (make step not run?)', id)); continue; end
            w = load(matf);
            if isfield(w,'out'), od = w.out; elseif isfield(w,'outDat'), od = w.outDat; else, error('no outDat'); end
            if isfield(od, 'moreThan1'), logmsg(logf, sprintf('SKIP_thresh %s (already final; make not run?)', id)); continue; end

            P = applyParams('threshTask', S.id);
            [outDat, raw, TTL] = assemble_outDat_all(S, P);
            outDat = downsample_data(outDat, P.fs_target);
            if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
            outDat = preprocess_macros(outDat, P);
            R = preprocess_respiration_wholetrace(outDat);

            % TASK-SPECIFIC: rebuild the 45-trial sniff-TTL table (start = sniff-1000)
            trial = 1:45;
            sniff = outDat.TTL.sniff;
            outDat.TTL = table(sniff(:)-1000, trial(:), sniff(:), 'variablenames', {'start','trial','sniff'});
            sniffs = detect_sniffs_from_TTLs(R, P, outDat);

            outDat.rspIDX = P.rspIDX; outDat.rspFlip = P.rspFlip;
            outDat.behDat = build_behavior_table_threshTask(sniffs, raw.beh);
            outDat = refine_onsets_with_phase(outDat, R, P);
            plot_sniff_epochs(outDat, R);
            outDat.moreThan1 = 0;

            if ~exist(preDir, 'dir'), mkdir(preDir); end
            save(fullfile(preDir, [S.id '_PEA_threshold_preproc.mat']), 'outDat', '-v7.3');
            hasMan = ismember('manOnset', outDat.behDat.Properties.VariableNames);
            logmsg(logf, sprintf('SUCCESS_thresh %s (manOnset=%d nBeh=%d fs=%d)', id, hasMan, height(outDat.behDat), outDat.fs));
            reproc_markDone(donePath, id);
        catch ME
            logmsg(logf, sprintf('FAIL_thresh %s : %s', id, ME.message));
        end
        clear outDat raw R sniffs TTL P od w
    end
    logmsg(logf, 'RUN_THRESH_MAIN_DONE'); fclose(logf);
end

function logmsg(fid, msg)
    fprintf('%s\n', msg);
    if fid > 0, fprintf(fid, '%s\n', msg); end
end
