function reproc_breathing_main(idxList)
% REPROC_BREATHING_MAIN  Portable main step (intermediate -> final) for
% breathingTask, mirroring breathingTaskPreProc_main.m (per-breath bmObj, block
% tagging, target-trace alignment, ECG/HRV, breath QC). Run AFTER
% reproc_backup('breathingTask') + reproc_breathing_make. Curated-only, resumes
% via done-marker, per-session failure isolation.
%
%   reproc_breathing_main()         % all curated breathing sessions
%   reproc_breathing_main(idxList)  % only these indices

    L = labPaths();
    addpath(genpath(L.repo));
    addpath(genpath(L.slowBreathing));
    addpath(genpath(L.eeglab));
    targTraceDir = L.targTraceDir;
    figPath = L.figPath;
    codePre = L.codePre;
    EEGLOC  = readtable(L.eegLocCsv);
    set(0, 'defaultfigurevisible', 'off');

    logf     = fopen(fullfile(reproc_root(), 'breathing_run.log'), 'a');
    donePath = fullfile(reproc_root(), 'breathing_done.txt');
    doneIds  = reproc_readDone(donePath);

    cfg = applyParams('breathingTask', 'main');
    if nargin < 1 || isempty(idxList), idxList = find(~strcmpi(cfg.paramSource, 'guess')); end
    logmsg(logf, sprintf('==== reproc_breathing_main start %s ====', datestr(now)));

    for s = idxList(:)'
        id = cfg.sessionIDs{s};
        if strcmpi(cfg.paramSource{s}, 'guess'), logmsg(logf, sprintf('GUESS_SKIP_breathing %s', id)); continue; end
        if any(strcmpi(id, doneIds)),            logmsg(logf, sprintf('ALREADY_DONE_breathing %s', id)); continue; end
        try
            S = struct('id', id, 'root', cfg.root{s}, 'fig', fullfile(figPath, id));
            preDir = fullfile(S.root, S.id, 'preProc');
            matf   = fullfile(preDir, [S.id '_breathingPreproc.mat']);
            if ~exist(matf, 'file'), logmsg(logf, sprintf('NOINT_breathing %s (make not run?)', id)); continue; end
            wv = load(matf); fn = fieldnames(wv); od = wv.(fn{1});
            if isfield(od, 'baseEmotion'), logmsg(logf, sprintf('SKIP_breathing %s (already final; make not run?)', id)); continue; end

            P   = applyParams('breathingTask', S.id);
            raw    = assembleRaw_breathingTask(S);
            outDat = assembleOutDat(raw, S, P);
            outDat = downsample_data(outDat, P.fs_target);
            if P.hasEEG,    outDat = preprocess_eeg(outDat, EEGLOC, P); end
            if P.hasMacros, outDat = preprocess_macros(outDat, P); end

            outDat = process_respiration_breathing(outDat, P);
            if isfield(outDat, 'TTL')
                TaskBreaks = [outDat.TTL/outDat.fs size(outDat.data,2)/outDat.fs];
            else
                TaskBreaks = 0:300:max(outDat.behDat.order)*300-10;
            end
            for cndi = 1:length(TaskBreaks)
                outDat.bmObj(outDat.bmObj(:,2) > TaskBreaks(cndi), 12) = cndi;
            end
            outDat.moreThan1 = 1;
            outDat.rspIDX = P.rspIDX; outDat.rspFlip = P.rspFlip;

            if ~strcmp(S.id, '250811_Dupi_NMH_TPB_1')
                [outDat, ~] = alignTargetBreathingTraceSimplify(outDat, targTraceDir);
            end
            outDat = build_behavior_table_breathingTask(outDat, outDat.bmObj);
            outDat = processECG(outDat, P);
            outDat = flagBadBreaths(outDat);

            R = preprocess_respiration_wholetrace(outDat);
            plot_sniff_epochs(outDat, R);
            plotBreathLengths(outDat, R);

            try
                csvDir = fullfile(codePre, 'closed-loop-respiration', 'processedBehavior');
                if isfolder(csvDir)
                    writetable(outDat.behDat, fullfile(csvDir, [outDat.sessID '_processedBreathing.csv']));
                end
            catch, end

            if ~exist(preDir, 'dir'), mkdir(preDir); end
            parSave(fullfile(preDir, [S.id '_breathingPreproc.mat']), outDat);
            hasMan = ismember('manOnset', outDat.behDat.Properties.VariableNames);
            logmsg(logf, sprintf('SUCCESS_breathing %s (manOnset=%d baseEmotion=%d nBeh=%d fs=%d)', ...
                id, hasMan, isfield(outDat,'baseEmotion'), height(outDat.behDat), outDat.fs));
            reproc_markDone(donePath, id);
        catch ME
            logmsg(logf, sprintf('FAIL_breathing %s : %s', id, ME.message));
        end
        clear outDat raw R wv od P
        close all force
    end
    logmsg(logf, 'RUN_BREATHING_MAIN_DONE'); fclose(logf);
end

function logmsg(fid, msg)
    fprintf('%s\n', msg);
    if fid > 0, fprintf(fid, '%s\n', msg); end
end
