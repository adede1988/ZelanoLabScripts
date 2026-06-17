function run_thresh_main(idxList)
% Robust runner mirroring threshPreProc_main.m (new loader + unchanged downstream).
% Processes only FRESH intermediates (no moreThan1). try/catch isolates failures.

    codePre = 'C:\Users\Adam\Documents\GitHub\';
    addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'));
    addpath(genpath([codePre 'ZelanoLabScripts']));
    figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';
    EEGLOC  = readtable(fullfile(codePre,'ZelanoLabScripts','eegLocs_standard_coords.csv'));
    set(0, 'defaultfigurevisible', 'off');

    cfg = applyParams('threshTask','main');
    if nargin < 1 || isempty(idxList), idxList = 1:numel(cfg.sessionIDs); end

    for s = idxList(:)'
        id = cfg.sessionIDs{s};
        try
            S = struct('id', id, 'root', cfg.root{s}, 'fig', fullfile(figPath, id));
            preDir = fullfile(S.root, S.id, 'preProc');
            matf = fullfile(preDir, [S.id '_PEA_threshold_preproc.mat']);
            if ~exist(matf, 'file'), fprintf('NOINT_thresh %s\n', id); continue; end
            w = load(matf);
            if isfield(w,'out'), od = w.out; elseif isfield(w,'outDat'), od = w.outDat; else, error('no outDat'); end
            if isfield(od,'moreThan1'), fprintf('SKIP_thresh %s (already final)\n', id); continue; end

            P = applyParams('threshTask', S.id);
            [outDat, raw, TTL] = assemble_outDat_all(S, P);
            outDat = downsample_data(outDat, P.fs_target);
            if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
            outDat = preprocess_macros(outDat, P);
            R = preprocess_respiration_wholetrace(outDat);

            trial = 1:45;
            sniff = outDat.TTL.sniff;
            x = table(sniff(:)-1000, trial(:), sniff(:), 'variablenames', {'start','trial','sniff'});
            outDat.TTL = x;
            sniffs = detect_sniffs_from_TTLs(R, P, outDat);

            outDat.rspIDX = P.rspIDX; outDat.rspFlip = P.rspFlip;
            outDat.behDat = build_behavior_table_threshTask(sniffs, raw.beh);
            outDat = refine_onsets_with_phase(outDat, R, P);
            plot_sniff_epochs(outDat, R);
            outDat.moreThan1 = 0;
            if ~exist(preDir,'dir'), mkdir(preDir); end
            save(fullfile(preDir, [S.id '_PEA_threshold_preproc.mat']), 'outDat', '-v7.3');
            fprintf('SUCCESS_thresh %s\n', id);
        catch ME
            fprintf('FAIL_thresh %s : %s\n', id, ME.message);
        end
        clear outDat raw R sniffs TTL P od w x   % free multi-GB session data before next iter
    end
    fprintf('RUN_THRESH_MAIN_DONE\n');
end
