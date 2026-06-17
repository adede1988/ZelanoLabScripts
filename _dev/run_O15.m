function run_O15(idxList, overwrite)
% Robust batch runner for the O15 pipeline. Mirrors O15PreProc_main.m exactly
% (same loader: applyParams + assemble_outDat_all; same downstream), but wraps
% each session in try/catch so one failure does not abort the batch, and logs a
% per-session SUCCESS/FAIL summary. overwrite=true reprocesses even if the final
% _O15preproc.mat already exists.

    if nargin < 2 || isempty(overwrite), overwrite = true; end

    codePre = 'C:\Users\Adam\Documents\GitHub\';
    addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'));
    addpath(genpath([codePre 'ZelanoLabScripts']));
    figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';
    EEGLOC  = readtable([codePre 'ZelanoLabScripts/eegLocs_standard_coords.csv']);
    set(0, 'defaultfigurevisible', 'off');

    cfg = applyParams('O15', 'main');
    if nargin < 1 || isempty(idxList), idxList = 1:numel(cfg.sessionIDs); end

    for s = idxList(:)'
        id = cfg.sessionIDs{s};
        try
            S = struct('id', id, 'root', cfg.root{s}, 'figPath', figPath);
            finalMat = fullfile(S.root, S.id, 'preProc', [S.id '_O15preproc.mat']);
            if ~overwrite && exist(finalMat, 'file')
                fprintf('SKIP %s (final exists)\n', id); continue;
            end
            if ~isfolder(fullfile(S.figPath, S.id)), mkdir(fullfile(S.figPath, S.id)); end

            P = applyParams('O15', S.id);
            [outDat, raw, TTL] = assemble_outDat_all(S, P);
            outDat.TTL = TTL;

            outDat = downsample_data(outDat, P.fs_target);
            if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
            outDat = preprocess_macros(outDat, P);

            outDat.moreThan1 = 1;
            outDat.rspIDX = P.rspIDX;
            outDat.rspFlip = P.rspFlip;

            R = preprocess_respiration_wholetrace(outDat);
            sniffs = detect_sniffs_from_TTLs(R, P, outDat);
            outDat.behDat = build_behavior_table_O15(sniffs, raw.beh);
            outDat = refine_onsets_with_phase(outDat, R, P);
            plot_sniff_epochs(outDat, R);

            if ~isfolder(fullfile(outDat.OGdataDir, 'preProc'))
                mkdir(fullfile(outDat.OGdataDir, 'preProc'));
            end
            save(fullfile(outDat.OGdataDir, 'preProc', [outDat.sessID '_O15preproc.mat']), ...
                 'outDat', '-v7.3');
            fprintf('SUCCESS_O15 %s\n', id);
        catch ME
            fprintf('FAIL_O15 %s : %s\n', id, ME.message);
        end
        clear outDat raw R sniffs TTL P   % free multi-GB session data before next iter
    end
    fprintf('RUN_O15_DONE\n');
end
