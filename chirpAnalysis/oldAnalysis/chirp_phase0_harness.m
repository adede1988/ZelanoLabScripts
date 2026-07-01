function results = chirp_phase0_harness(sweep, C)
% CHIRP_PHASE0_HARNESS  Validation sweep (spec 5) -- GATES interpretation of the real-data battery.
%   results = chirp_phase0_harness(sweep, C)
%   sweep fields (all optional; sensible defaults):
%     .overlapGrid [ms]  .DfGrid [Hz]  .ampGrid {[hi lo]...}  .families {..}
%     .nTrials  .fLo  .snrDb  .burstMs  .fastDropMs  .windowSources {'pipeline','oracle'}
%     .outDir  (default <groupDir>/chirpAnalysis/phase0)
%
%   Synthesizes trials, pushes them through the REAL battery (chirp_tfr_faslt -> chirp_ridge ->
%   phase/beat/chirplet/temporal) on BOTH ground-truth ('oracle') and ridge-derived ('pipeline')
%   windows. Emits a tidy per-cell table (phase0_cells.csv) + results .mat. The asymmetric reads:
%   phase concentrated -> single; beat significant -> dual; chirplet fracSingle/fracDual; temporal
%   corr -> coupling. Deliverable curves (beat-min-overlap, chirplet single-bias) are derived from
%   the table downstream.

    if nargin < 2 || isempty(C), C = chirp_config(); end
    d = defaults();
    if nargin >= 1 && ~isempty(sweep)
        fn = fieldnames(sweep); for i=1:numel(fn), d.(fn{i}) = sweep.(fn{i}); end
    end
    outDir = d.outDir; if isempty(outDir), outDir = fullfile(resolveOutDir(C), 'phase0'); end
    if ~isfolder(outDir), mkdir(outDir); end

    fs = C.fs;
    [off, coreIdx, tMs] = padGrid(C);

    rows = {};
    fprintf('Phase-0 harness: %d families x %d overlaps x %d Df x %d amp x %d windowSrc\n', ...
        numel(d.families), numel(d.overlapGrid), numel(d.DfGrid), numel(d.ampGrid), numel(d.windowSources));
    for fam = d.families
        family = fam{1};
        isSingle = startsWith(family,'single');
        ampList = d.ampGrid; if isSingle, ampList = {NaN}; end   % single: amp irrelevant
        ovList = d.overlapGrid; if isSingle, ovList = unique(min(d.burstMs, d.overlapGrid)); end
        for Df = d.DfGrid
            for ai = 1:numel(ampList)
                amp = ampList{ai};
                for ov = ovList
                    % synth a cell of trials
                    X = zeros(d.nTrials, numel(off)); GT = cell(1,d.nTrials);
                    for k = 1:d.nTrials
                        P = struct('overlapMs',ov,'Df',Df,'ampRatio',amp,'fLo',d.fLo, ...
                            'snrDb',d.snrDb,'burstMs',d.burstMs,'fastDropMs',d.fastDropMs, ...
                            'seed', cellSeed(C.rngSeed, family, Df, ai, ov, k));
                        [X(k,:), GT{k}] = chirp_synth_trial(family, P, C);
                    end
                    E = struct('dataPad',X,'coreIdx',coreIdx,'tMs',tMs,'valid',true(d.nTrials,1), ...
                        'fs',fs,'onsets',(1:d.nTrials)','padSec',C.padSec,'epochWin',C.epochWin);
                    tfr = chirp_tfr_faslt(E, C);
                    R0  = chirp_ridge(tfr, E, C);
                    for ws = d.windowSources
                        Rw = R0;
                        if strcmp(ws{1},'oracle')
                            for k = 1:d.nTrials
                                Rw.trial(k).transitionWin = GT{k}.transitionWinMs;
                                Rw.trial(k).coexistWin    = GT{k}.coexistWinMs;
                            end
                        end
                        m = runBattery(E, Rw, tfr, tMs, C);
                        rows(end+1,:) = {family, isSingle, Df, ampStr(amp), ov, ws{1}, ...
                            m.perm_z, m.perm_p, m.Rbar, m.perm_z0, m.vtest_p, m.meanResidDeg, m.crossR, m.crossP, m.phaseNUsed, ...
                            m.beat_fracSig, m.beat_medLogSNR, m.beat_nUsed, ...
                            m.fracSingle, m.fracDual, m.fracAmbig, ...
                            m.corrHiLo, m.varLatDiff, m.meanNRidges, m.fracNR2}; %#ok<AGROW>
                    end
                    fprintf('  %-10s Df=%2d amp=%-4s ov=%4d done\n', family, Df, ampStr(amp), ov);
                end
            end
        end
    end

    T = cell2table(rows, 'VariableNames', {'family','isSingle','Df','ampRatio','overlapMs', ...
        'windowSource','perm_z','perm_p','Rbar','perm_z0','vtest_p','meanResidDeg','crossR','crossP','phaseNUsed','beat_fracSig','beat_medLogSNR', ...
        'beat_nUsed','frac_single','frac_dual','frac_ambig','corrHiLo','varLatDiff','meanNRidges','fracNR2'});
    writetable(T, fullfile(outDir,'phase0_cells.csv'));
    results = struct('T',T,'sweep',d,'C',C);
    save(fullfile(outDir,'phase0_results.mat'),'results','-v7.3');
    fprintf('Phase-0 done -> %s\n', fullfile(outDir,'phase0_cells.csv'));
    try, phase0_plots(T, outDir); catch ME, fprintf('plot skipped: %s\n', ME.message); end
end

% ===== battery on one cell =====
function m = runBattery(E, R, tfr, tMs, C)
    ph = chirp_phase_test(E, R, C);
    be = chirp_beat_test(E, R, C);
    td = chirp_temporal_decoup(E, R, C);
    nTr = size(E.dataPad,1); cls = cell(1,nTr);
    for i = 1:nTr
        wfit = R.trial(i).transitionWin; if isempty(wfit), wfit = R.trial(i).coexistWin; end
        o = chirp_chirplet(E.dataPad(i,E.coreIdx), tMs, wfit, C);
        cls{i} = o.classification;
    end
    nrs = arrayfun(@(t) t.nRidgesSig, R.trial);
    m.perm_z = ph.summary.perm_z; m.perm_p = ph.summary.perm_p; m.Rbar = ph.summary.Rbar;
    m.perm_z0 = ph.summary.perm_z0; m.vtest_p = ph.summary.vtest_p_vs0; m.meanResidDeg = ph.summary.meanResidDeg;
    m.crossR = ph.summary.crossTrialR; m.crossP = ph.summary.crossTrial_p;
    m.phaseNUsed = ph.summary.nUsed;
    m.beat_fracSig = ph_nan(be.summary.fracSig); m.beat_medLogSNR = be.summary.medianLogSNR;
    m.beat_nUsed = be.summary.nUsed;
    m.fracSingle = mean(strcmp(cls,'single')); m.fracDual = mean(strcmp(cls,'dual'));
    m.fracAmbig = mean(strcmp(cls,'ambig'));
    m.corrHiLo = td.summary.corrHiLo; m.varLatDiff = td.summary.varLatDiff;
    m.meanNRidges = mean(nrs,'omitnan'); m.fracNR2 = mean(nrs>1);
end

% ===== helpers =====
function d = defaults()
    d.overlapGrid = [50 100 200 500 1000 2000];
    d.DfGrid = [3 5 7 10];
    d.ampGrid = {[1 1],[2 1],[4 1]};
    d.families = {'singleLin','singleNL','dualIndep','dualLocked'};
    d.nTrials = 30; d.fLo = 33; d.snrDb = 8; d.burstMs = 1500; d.fastDropMs = 400;
    d.windowSources = {'pipeline','oracle'}; d.outDir = '';
end
function [off, coreIdx, tMs] = padGrid(C)
    fs = C.fs;
    s0p = round((C.epochWin(1)-C.padSec)*fs); s1p = round((C.epochWin(2)+C.padSec)*fs);
    off = s0p:s1p; coreIdx = off>=round(C.epochWin(1)*fs) & off<=round(C.epochWin(2)*fs);
    tMs = off(coreIdx)/fs*1000;
end
function s = cellSeed(base, family, Df, ai, ov, k)
    s = mod(base + 1000*sum(double(family)) + 137*Df + 31*ai + 7*ov + k, 2^31-1);
end
function a = ampStr(amp), if isscalar(amp)&&isnan(amp), a='NA'; else, a=sprintf('%d:%d',amp(1),amp(2)); end, end
function v = ph_nan(x), if isempty(x), v=NaN; else, v=x; end, end
function o = resolveOutDir(C)
    if ~isempty(C.outDir), o = C.outDir; return; end
    try, L = labPaths(); o = fullfile(L.figPath,'groupStatFigs','chirpAnalysis'); catch, o = pwd; end
end

function phase0_plots(T, outDir)
% beat-minimum-overlap (dualIndep, pipeline) + chirplet single-bias (dualIndep) vs overlap, per Df
    pdat = T(strcmp(T.family,'dualIndep') & strcmp(T.windowSource,'pipeline'), :);
    if isempty(pdat), return; end
    Dfs = unique(pdat.Df);
    fh = figure('Visible','off','Position',[80 80 1000 420]);
    subplot(1,2,1); hold on;
    for q = 1:numel(Dfs)
        s = pdat(pdat.Df==Dfs(q),:); s = sortrows(s,'overlapMs');
        plot(s.overlapMs, s.beat_fracSig, '-o', 'DisplayName', sprintf('Df=%d',Dfs(q)));
    end
    xlabel('overlap (ms)'); ylabel('beat fracSig'); title('Beat sensitivity vs overlap (dualIndep)');
    legend('Location','best'); grid on; ylim([0 1]);
    subplot(1,2,2); hold on;
    for q = 1:numel(Dfs)
        s = pdat(pdat.Df==Dfs(q),:); s = sortrows(s,'overlapMs');
        plot(s.overlapMs, s.frac_single, '-o', 'DisplayName', sprintf('Df=%d',Dfs(q)));
    end
    xlabel('overlap (ms)'); ylabel('chirplet frac SINGLE (on true dual)');
    title('Chirplet single-bias vs overlap (should drop)'); legend('Location','best'); grid on; ylim([0 1]);
    exportgraphics(fh, fullfile(outDir,'phase0_curves.png'), 'Resolution',120); close(fh);
end
