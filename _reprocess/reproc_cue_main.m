function reproc_cue_main(idxList)
% REPROC_CUE_MAIN  Portable main step (intermediate -> final) for cueTask,
% mirroring cueTaskPreProc_main.m. Run AFTER reproc_backup('cueTask') +
% reproc_cue_make (which regenerate the intermediates from raw). Processes only
% CURATED sessions, resumes via a done-marker, isolates per-session failures
% (incl. ZLP:blinkAmbiguous -> skipped + review plot saved by blinkRemoveWrapper).
%
%   reproc_cue_main()         % all curated cue sessions
%   reproc_cue_main(idxList)  % only these indices

    L = labPaths();
    addpath(genpath(L.repo));
    addpath(genpath(L.eeglab));
    figPath = L.figPath;
    EEGLOC  = readtable(L.eegLocCsv);
    set(0, 'defaultfigurevisible', 'off');

    logf     = fopen(fullfile(reproc_root(), 'cue_run.log'), 'a');
    donePath = fullfile(reproc_root(), 'cue_done.txt');
    doneIds  = reproc_readDone(donePath);

    cfg = applyParams('cueTask', 'main');
    if nargin < 1 || isempty(idxList), idxList = find(~strcmpi(cfg.paramSource, 'guess')); end
    logmsg(logf, sprintf('==== reproc_cue_main start %s ====', datestr(now)));

    for s = idxList(:)'
        id = cfg.sessionIDs{s};
        if strcmpi(cfg.paramSource{s}, 'guess'), logmsg(logf, sprintf('GUESS_SKIP_cue %s', id)); continue; end
        if any(strcmpi(id, doneIds)),            logmsg(logf, sprintf('ALREADY_DONE_cue %s', id)); continue; end
        try
            S = struct('id', id, 'root', cfg.root{s}, 'fig', fullfile(figPath, id));
            preDir = fullfile(S.root, S.id, 'preProc');
            matf   = fullfile(preDir, [S.id '_cueTaskPreProc.mat']);
            if ~exist(matf, 'file'), logmsg(logf, sprintf('NOINT_cue %s (make step not run?)', id)); continue; end
            w = load(matf);
            if isfield(w,'out'), od = w.out; elseif isfield(w,'outDat'), od = w.outDat; else, error('no outDat'); end
            if isfield(od, 'moreThan1'), logmsg(logf, sprintf('SKIP_cue %s (already final; make not run?)', id)); continue; end

            P = applyParams('cueTask', S.id);
            [outDat, raw, TTL] = assemble_outDat_all(S, P);
            outDat = downsample_data(outDat, P.fs_target);
            if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
            outDat = preprocess_macros(outDat, P);
            R = preprocess_respiration_wholetrace(outDat);
            sniffs = detect_sniffs_from_TTLs(R, P, outDat);
            outDat.moreThan1 = 0;
            outDat.rspIDX = P.rspIDX; outDat.rspFlip = P.rspFlip;
            outDat.behDat = build_behavior_table_cueTask(sniffs, raw.beh);
            outDat = refine_onsets_with_phase(outDat, R, P);
            plot_sniff_epochs(outDat, R);

            if ~exist(preDir, 'dir'), mkdir(preDir); end
            save(fullfile(preDir, [S.id '_cueTaskPreproc.mat']), 'outDat', '-v7.3');
            hasMan = ismember('manOnset', outDat.behDat.Properties.VariableNames);
            logmsg(logf, sprintf('SUCCESS_cue %s (manOnset=%d nBeh=%d fs=%d)', id, hasMan, height(outDat.behDat), outDat.fs));
            reproc_markDone(donePath, id);
        catch ME
            logmsg(logf, sprintf('FAIL_cue %s : %s', id, ME.message));
        end
        clear outDat raw R sniffs TTL P od w
    end
    logmsg(logf, 'RUN_CUE_MAIN_DONE'); fclose(logf);
end

function logmsg(fid, msg)
    fprintf('%s\n', msg);
    if fid > 0, fprintf(fid, '%s\n', msg); end
end
