function run_chirp_analysis_pipeline(sessFilter, opts)
% RUN_CHIRP_ANALYSIS_PIPELINE  Phase-1 driver: per-session -> per-channel -> per-trial battery.
%   run_chirp_analysis_pipeline()                 % all FRESH cue sessions on E:
%   run_chirp_analysis_pipeline({'<sessID>'})     % one/few sessions
%   run_chirp_analysis_pipeline([], struct('profile',true,'maxTrials',8,'saveFigs',true))
%
%   For each fresh cue final (chirp_session_table, D5): load bestMac (OB) + all macBP; epoch on
%   behDat.finalOnset ([-1,+3]s, padded); FASLT TFR -> ridge -> phase/beat/chirplet/temporal
%   (bestMac) + spatial (across macBP). Assembles the spec-7.1 sub struct (saved as
%   <sessID>_chirp.mat) and appends the spec-7.2 tidy per-trial CSVs. Figures under
%   <figRoot>/<id>/singleTrialTF/. Outputs to E: via CHIRP_OUTDIR/CHIRP_FIGROOT.
%
%   opts: .profile(false) .maxTrials(inf) .saveFigs(false in batch) .nSurr .ridgeNSurr .useParfor

    if nargin<1, sessFilter = []; end
    if nargin<2, opts = struct(); end
    setup_chirpAnalysis_paths();
    ov = struct(); if isfield(opts,'nSurr'), ov.surr.n=opts.nSurr; end
    if isfield(opts,'ridgeNSurr'), ov.ridge.nSurr=opts.ridgeNSurr; end
    C = chirp_config(ov);

    profile = getf(opts,'profile',false);
    maxTrials = getf(opts,'maxTrials',inf);
    saveFigs = getf(opts,'saveFigs', profile);

    outDir = resolveOutDir(C); if ~isfolder(outDir), mkdir(outDir); end
    figRoot = C.figRoot; if isempty(figRoot), try, L=labPaths(); figRoot=L.figPath; catch, figRoot=outDir; end, end
    fprintf('chirp pipeline: outDir=%s\n figRoot=%s\n', outDir, figRoot);

    freshOnly = ~getf(opts,'allowStale',false);       % D5: fresh by default; allowStale=mechanics test
    T = chirp_session_table(C, freshOnly, true);
    if ~freshOnly, fprintf('  *** allowStale=true: STALE finals permitted -- MECHANICS TEST ONLY, not a result ***\n'); end
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end
    if isfield(opts,'maxSessions') && height(T) > opts.maxSessions, T = T(1:opts.maxSessions, :); end
    if isempty(T), fprintf('No sessions to process (no fresh finals on E: yet -- sweep later).\n'); return; end

    summaryRows = {};
    for si = 1:height(T)
        id = char(T.sessID(si)); fp = char(T.path(si)); t0 = tic;
        fprintf('\n==== %d/%d  %s ====\n', si, height(T), id);
      try   % per-session guard: a single bad session must not abort the whole batch
        S = chirp_load_session(fp, C);
        if ~S.ok, fprintf('  SKIP: %s\n', S.msg); continue; end
        if ~isempty(S.msg), fprintf('  note: %s\n', S.msg); end

        E   = chirp_epoch(S.bestSig, S.fs, S.onsets, C);
        tfr = chirp_tfr_faslt(E, C);
        R   = chirp_ridge(tfr, E, C);
        ph  = chirp_phase_test(E, R, C);
        be  = chirp_beat_test(E, R, C);
        td  = chirp_temporal_decoup(E, R, C);
        nTr = size(E.dataPad,1);
        chOut = repmat(chirp_chirplet([],[],[],C), nTr, 1);   % blank
        for i = 1:nTr
            if i>maxTrials, break; end
            wfit = R.trial(i).transitionWin; if isempty(wfit), wfit = R.trial(i).coexistWin; end
            chOut(i) = chirp_chirplet(E.dataPad(i,E.coreIdx), E.tMs, wfit, C);
        end
        sp = chirp_spatial(S.macSig, S.fs, S.onsets, S.goodMac, R, C);
        bl = chirp_measure_burst(S, C);

        sub = assembleSub(id, T(si,:), S, C, R, ph, be, td, chOut, sp, bl);
        save(fullfile(outDir, [id '_chirp.mat']), 'sub', '-v7.3');

        appendTrialCSVs(outDir, id, T(si,:), S, R, ph, be, td, chOut, sp);
        summaryRows(end+1,:) = summaryRow(id, T(si,:), S, ph, be, td, chOut, sp, bl, R); %#ok<AGROW>

        if saveFigs
            figDir = fullfile(figRoot, id, C.figSub);
            chirp_fig_montage(tfr, R, S.bestMac, id, S.bestIdx, figDir, C);
            np = min(nTr, getf(opts,'nFigTrials',12));
            for i = round(linspace(1,nTr,np))
                if i<=numel(chOut), chirp_fig_trial(tfr, R, S.bestMac, i, id, S.bestIdx, figDir, C, chOut(i)); end
            end
        end
        fprintf('  done in %.1f s | phase perm_z=%.2f | beat fracSig=%.2f | chirp S/D=%.2f/%.2f | nRidge2=%.2f\n', ...
            toc(t0), ph.summary.perm_z, naN(be.summary.fracSig), ...
            mean(strcmp({chOut.classification},'single')), mean(strcmp({chOut.classification},'dual')), mean(arrayfun(@(t)t.nRidgesSig,R.trial)>1));
      catch ME
        fprintf('  *** SESSION FAILED (%s): %s\n', id, ME.message);
        fprintf('      %s\n', getReport(ME, 'basic', 'hyperlinks','off'));
      end
        clear S E tfr R ph be td chOut sp bl sub;
    end

    if ~isempty(summaryRows)
        Sm = cell2table(summaryRows, 'VariableNames', summaryVars());
        writetable(Sm, fullfile(outDir, 'subject_summary.csv'));
        fprintf('\nWrote subject_summary.csv (%d sessions) -> %s\n', height(Sm), outDir);
    end
end

% ================= assembly =================
function sub = assembleSub(id, Trow, S, C, R, ph, be, td, chOut, sp, bl)
    sub = struct();
    sub.sessID = id; sub.group = char(Trow.group); sub.type = char(Trow.type);
    sub.bestMac = S.bestMac; sub.macLabs = {S.macLabs}; sub.fs = S.fs; sub.cfg = C;
    sub.ridge   = struct('trial', R.trial, 'anchors', R.anchors);
    sub.phaseTest = ph;
    sub.beating   = be;
    sub.temporalDecoup = td;
    sub.chirplet = struct('trial', chOut, 'summary', chirpletSummary(chOut));
    sub.spatial  = sp;
    sub.burst    = bl.B; sub.snrDb = bl.snrDb;
end

function s = chirpletSummary(chOut)
    cls = {chOut.classification};
    c1 = arrayfun(@(o) o.c1, chOut);
    traj = arrayfun(@(o) o.trajConnectivity, chOut);
    s = struct('medianC1', median(c1,'omitnan'), 'fracDownchirp', mean(c1<0), ...
        'fracSingleTraj', mean(strcmp(cls,'single')), 'fracDualTraj', mean(strcmp(cls,'dual')), ...
        'fracAmbig', mean(strcmp(cls,'ambig')), 'medianTraj', median(traj,'omitnan'));
end

% ================= CSV append =================
function appendTrialCSVs(outDir, id, Trow, S, R, ph, be, td, chOut, sp)
    grp = char(Trow.group); typ = char(Trow.type);
    nTr = numel(R.trial);
    key = @(i) {string(id), string(grp), string(typ), string(S.bestMac), i};
    % ridge
    rr = {};
    for i=1:nTr, t=R.trial(i);
        rr(end+1,:) = [key(i), {t.nRidgesSig, t.ridgeOverlap, t.powerDipStat, ...
            winv(t.transitionWin,1), winv(t.transitionWin,2), winv(t.coexistWin,1), winv(t.coexistWin,2), t.included}]; %#ok<AGROW>
    end
    appendCsv(fullfile(outDir,'ridge_trials.csv'), rr, ...
        {'session','group','type','channel','trial','nRidgesSig','ridgeOverlap','powerDipStat','transStartMs','transEndMs','coexistStartMs','coexistEndMs','included'});
    % phase
    pr = {};
    for i=1:nTr, t=ph.trial(i);
        pr(end+1,:) = [key(i), {nanmeanv([t.meanRbar]), numel(t.subWinResid), numel(t.slipLatencies), t.included}]; %#ok<AGROW>
    end
    appendCsv(fullfile(outDir,'phaseTest_trials.csv'), pr, ...
        {'session','group','type','channel','trial','meanRbar','nSubWin','nSlips','included'});
    % beat
    br = {};
    for i=1:nTr, t=be.trial(i);
        br(end+1,:) = [key(i), {t.fBeat, t.modSNR, t.modIndex, t.nBeatCycles, t.deepNullPeriodicity, t.modSNR_baseline, t.modSNR_matchedPower, t.sigSurr, t.included}]; %#ok<AGROW>
    end
    appendCsv(fullfile(outDir,'beating_trials.csv'), br, ...
        {'session','group','type','channel','trial','fBeat','modSNR','modIndex','nBeatCycles','deepNull','modSNR_baseline','modSNR_matched','sigSurr','included'});
    % chirplet
    cr = {};
    for i=1:nTr, o=chOut(i);
        cr(end+1,:) = [key(i), {o.nAtomsSig, o.trajConnectivity, string(o.classification), o.c1, o.dBIC, o.f1, o.f2, o.tSep, o.overlap, o.included}]; %#ok<AGROW>
    end
    appendCsv(fullfile(outDir,'chirplet_trials.csv'), cr, ...
        {'session','group','type','channel','trial','nAtomsSig','trajConnectivity','classification','c1','dBIC','f1','f2','tSep','overlap','included'});
    % temporal
    tr = {};
    for i=1:nTr, t=td.trial(i);
        tr(end+1,:) = [key(i), {t.latHi, t.latLo, t.latDiff, t.included}]; %#ok<AGROW>
    end
    appendCsv(fullfile(outDir,'temporalDecoup_trials.csv'), tr, ...
        {'session','group','type','channel','trial','latHi','latLo','latDiff','included'});
    % spatial (per trial; profileSim only -- profiles are vectors)
    sr = {};
    for i=1:numel(sp.trial), t=sp.trial(i);
        sr(end+1,:) = {string(id), string(grp), string(typ), i, sp.summary.eligible, sp.summary.nGoodContacts, t.profileSim, t.included}; %#ok<AGROW>
    end
    appendCsv(fullfile(outDir,'spatial_trials.csv'), sr, ...
        {'session','group','type','trial','eligible','nGoodContacts','profileSim','included'});
end

function appendCsv(fn, rows, vars)
    if isempty(rows), return; end
    Tt = cell2table(rows, 'VariableNames', vars);
    if isfile(fn), writetable(Tt, fn, 'WriteMode','append'); else, writetable(Tt, fn); end
end

% ================= subject summary =================
function r = summaryRow(id, Trow, S, ph, be, td, chOut, sp, bl, R)
    cs = arrayfun(@(o) o.classification, chOut, 'uni', 0);
    nr2 = mean(arrayfun(@(t) t.nRidgesSig, R.trial) > 1);
    r = {string(id), string(Trow.group), string(Trow.type), string(S.bestMac), numel(S.onsets), ...
        ph.summary.nUsed, ph.summary.Rbar, ph.summary.perm_z, ph.summary.perm_p, ph.summary.perm_z0, ph.summary.perm_p0, ph.summary.vtest_p_vs0, ph.summary.crossTrialR, ph.summary.crossTrial_p, ...
        be.summary.nUsed, be.summary.medianLogSNR, naN(be.summary.fracSig), ...
        median(arrayfun(@(o)o.c1,chOut),'omitnan'), mean(strcmp(cs,'single')), mean(strcmp(cs,'dual')), mean(strcmp(cs,'ambig')), ...
        td.summary.corrHiLo, td.summary.varLatDiff, td.summary.nUsed, ...
        sp.summary.eligible, sp.summary.nGoodContacts, sp.summary.meanProfileSim, sp.summary.nUsed, ...
        nr2, bl.snrDb, bl.B.summary.median, bl.B.summary.p90, bl.B.coverage};
end
function v = summaryVars()
    v = {'session','group','type','bestMac','nTrials', ...
        'phase_nUsed','phase_Rbar','phase_permZ','phase_permP','phase_permZ0','phase_permP0','phase_vtestP','phase_crossR','phase_crossP', ...
        'beat_nUsed','beat_medLogSNR','beat_fracSig', ...
        'chirp_medianC1','chirp_fracSingle','chirp_fracDual','chirp_fracAmbig', ...
        'temp_corrHiLo','temp_varLatDiff','temp_nUsed', ...
        'spatial_eligible','spatial_nGood','spatial_meanSim','spatial_nUsed', ...
        'ridge_fracNR2','snrDb','burst_medianMs','burst_p90Ms','burst_coverage'};
end

% ================= small helpers =================
function v = getf(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function v = naN(x), if isempty(x), v=NaN; else, v=x; end, end
function v = nanmeanv(x), if isempty(x), v=NaN; else, v=mean(x,'omitnan'); end, end
function v = winv(w,k), if isempty(w)||numel(w)<k, v=NaN; else, v=w(k); end, end
function o = resolveOutDir(C)
    if ~isempty(C.outDir), o=C.outDir; return; end
    try, L=labPaths(); o=fullfile(L.figPath,'groupStatFigs','chirpAnalysis'); catch, o=fullfile(pwd,'chirpOut'); end
end
