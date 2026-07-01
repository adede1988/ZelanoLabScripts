function chirpAnalysisV2(sessFilter, opts)
% CHIRPANALYSISV2  Driver for the V2 OB sniff-locked gamma ridge + phase-continuity analysis.
%   chirpAnalysisV2()                       % all fresh cue finals under C.dataRoot
%   chirpAnalysisV2({'<sessID>'})           % one/few sessions (profiling)
%   chirpAnalysisV2([], struct('saveFig',true,'saveMat',true,'maxSessions',1))
%
%   Per session (cueTask, spec chirpAnalysis_2.md):
%     load final -> bestMac + behDat.finalOnset -> cue_noise_trials QC (>=80% retained) ->
%     v2_tfr (epoch [-1,+3]s pad 1.5s, FASLT 20-70, baseline z-score) ->
%     v2_ridges (tfridge primary + Gaussian-FWHM peel -> secondary) -> ridgeInfo ->
%     v2_powerphase -> powerPhaseContinuity -> figures -> write both fields into outDat and
%     re-save the final (E:). Copyback E:->R: (finals + figs) is a separate step from home.
%
%   ALL helper functions are called from here (spec 2.6). Outputs on E:; figures under
%   C.figRootE\<sessID>\<figSub>. Everything task-agnostic except sniff selection (cue = cued).

    if nargin<1, sessFilter = []; end
    if nargin<2, opts = struct(); end
    setup_chirpAnalysis_paths(false);
    here = fileparts(mfilename('fullpath')); repo = fileparts(here);
    addpath(fullfile(repo,'cueAnalysis'));   % cue_noise_trials.m (spec 2.4 QC)
    C = v2_config();
    saveFig = getf(opts,'saveFig',true);
    saveMat = getf(opts,'saveMat',true);

    T0 = tic;
    files = listCueFinals(C);
    if isempty(files), fprintf('No cue finals under %s\n', C.dataRoot); return; end
    if ~isempty(sessFilter)
        keep = ismember({files.sessID}, cellstr(string(sessFilter)));
        files = files(keep);
    end
    if isfield(opts,'maxSessions') && numel(files) > opts.maxSessions
        files = files(1:opts.maxSessions);
    end
    fprintf('chirpAnalysisV2: %d session(s); figRootE=%s\n', numel(files), C.figRootE);

    outDir = getenvOr('CHIRP_V2OUT','E:\chirpV2out'); if ~isfolder(outDir), mkdir(outDir); end
    rows = {};

    for si = 1:numel(files)
        id = files(si).sessID; fp = files(si).path; t0 = tic;
        fprintf('\n==== %d/%d  %s ====\n', si, numel(files), id);
        try
            % ---- load full final (preserve top-level var name for save-back) ----
            s = load(fp); fn = fieldnames(s); vname = fn{1}; od = s.(vname); clear s;
            [bestSig, bestLab, onsets, behDat, fs, figR] = extractSubstrate(od, id);
            N = numel(onsets);

            % ---- noise QC on bestMac (spec 2.4) ----
            nq = cue_noise_trials(bestSig, fs, onsets, C.noise.epWin, C.noise.K, C.noise.winMs);
            good = nq.ok & ~nq.noisy;
            retention = mean(good);
            fprintf('  bestMac=%s  N=%d  retained=%.0f%% (K=%g)\n', bestLab, N, 100*retention, C.noise.K);

            % ---- TFR (epoch, FASLT, baseline z) ----
            T = v2_tfr(bestSig, fs, onsets, good, C);
            fprintf('  TFR done: %d freqs x %d samp x %d good trials\n', numel(T.freqs), T.nCore, nnz(T.good));

            % ---- ridges (primary + peel + secondary) ----
            ridgeInfo = v2_ridges(T, C);
            ridgeInfo.bestMac = bestLab; ridgeInfo.retention = retention; ridgeInfo.noiseK = C.noise.K;
            ridgeInfo.noisyMask = nq.noisy(:); ridgeInfo.okMask = nq.ok(:);

            % ---- power / phase continuity (+ ridgeBurst) ----
            pp = v2_powerphase(T, ridgeInfo, C);

            % ---- peak-aligned analysis (peakLockedRidge + freq/time-centered mean TFR) ----
            [plr, meanPeakTFR, relFreq, relTimeMs] = v2_peaklocked(T, pp, C);
            ridgeInfo.peakLockedRidge = plr;

            % ---- write new fields into outDat ----
            od.ridgeInfo = ridgeInfo;
            od.powerPhaseContinuity = pp;

            % ---- group + subject mean matrices for the group step ----
            grp = 'Dupi'; if isfield(od,'type')&&~isempty(od.type), grp = char(string(od.type)); end
            gidx = find(T.good);
            meanTFR = mean(T.zTFR(:,:,gidx), 3, 'omitnan');           % onset-locked subject mean
            [Mpk, phaseRelT] = v2_peakcenter(pp.peakPhaseConsistency,  pp.tMs, pp.gammaPeakTime, 1000, fs);
            [Mrg, ~        ] = v2_peakcenter(pp.ridgePhaseConsistency, pp.tMs, pp.gammaPeakTime, 1000, fs);
            peakConsMean  = mean(abs(Mpk(gidx,:)),1,'omitnan');
            ridgeConsMean = mean(abs(Mrg(gidx,:)),1,'omitnan');

            % ---- figures (E: mirror of the subject figs path) ----
            figDirE = fullfile(C.figRootE, id, C.figSub);
            if saveFig
                v2_figs_ridge(T, ridgeInfo, figDirE, id, bestLab, C);
                v2_figs_phase(pp, figDirE, id, C);
                v2_fig_peaktfr(meanPeakTFR, relFreq, relTimeMs, ...
                    fullfile(figDirE, sprintf('peakAlignedTFR_%s.png', id)), ...
                    sprintf('%s  peak-aligned mean TFR', strrep(id,'_','\_')));
            end

            % ---- per-subject aggregation (.mat) for the group step ----
            agg = struct('sessID',id,'group',grp,'bestMac',bestLab,'retention',retention, ...
                'meanTFR',meanTFR,'freqs',T.freqs,'tMs',T.tMs, ...
                'meanPeakTFR',meanPeakTFR,'relFreq',relFreq,'relTimeMs',relTimeMs, ...
                'phaseRelT',phaseRelT,'peakConsMean',peakConsMean,'ridgeConsMean',ridgeConsMean, ...
                'peakBurst',pp.peakBurstLength(~isnan(pp.peakBurstLength)), ...
                'ridgeBurst',pp.ridgeBurstLength(~isnan(pp.ridgeBurstLength)), ...
                'medPeakBurst',nanmed(pp.peakBurstLength),'medRidgeBurst',nanmed(pp.ridgeBurstLength));
            aggDir = fullfile(outDir,'agg'); if ~isfolder(aggDir), mkdir(aggDir); end
            save(fullfile(aggDir,[id '_agg.mat']),'agg','-v7.3');

            % ---- re-save the final on E: (overwrite; same top-level var) ----
            if saveMat
                tmp = struct(); tmp.(vname) = od; save(fp, '-struct', 'tmp', '-v7.3');
            end

            rows(end+1,:) = summaryRow(id, bestLab, N, retention, pp); %#ok<AGROW>
            fprintf('  done in %.1f s | medPeakBurst=%.0f medRidgeBurst=%.0f medPeakF=%.1fHz hasPeak=%.0f%%\n', ...
                toc(t0), nanmed(pp.peakBurstLength), nanmed(pp.ridgeBurstLength), nanmed(pp.gammaPeakFrequency), 100*nanmean(pp.hasPeak));
            clear T ridgeInfo pp od meanTFR meanPeakTFR Mpk Mrg;
        catch ME
            fprintf('  *** SESSION FAILED (%s): %s\n', id, ME.message);
            fprintf('      %s\n', getReport(ME,'basic','hyperlinks','off'));
        end
    end

    if ~isempty(rows)
        Sm = cell2table(rows, 'VariableNames', summaryVars());
        writetable(Sm, fullfile(outDir,'v2_subject_summary.csv'));
        fprintf('\nWrote v2_subject_summary.csv (%d sessions) -> %s\n', height(Sm), outDir);
    end

    % ---- group-level figures + R-stats CSVs (spec 4/5) ----
    if getf(opts,'doGroup', numel(files)>1)
        try
            v2_group(fullfile(outDir,'agg'), fullfile(C.figRootE,'groupStatFigs'), outDir, C);
        catch GE
            fprintf('group step failed: %s\n%s\n', GE.message, getReport(GE,'basic','hyperlinks','off'));
        end
    end
    fprintf('chirpAnalysisV2 total %.1f min\n', toc(T0)/60);
end

% =================== helpers (all local, called only from here) ===================
function files = listCueFinals(C)
    d = dir(fullfile(C.dataRoot,'**','*_cueTaskpreproc.mat'));   % case-insensitive on Windows
    files = struct('sessID',{},'path',{});
    for k = 1:numel(d)
        fp = fullfile(d(k).folder, d(k).name);
        if datetime(d(k).datenum,'ConvertFrom','datenum') < C.freshCutoff, continue; end
        [~,nm] = fileparts(d(k).name);
        id = regexprep(nm, '_cuetaskpreproc$', '', 'ignorecase');
        files(end+1) = struct('sessID',id,'path',fp); %#ok<AGROW>
    end
end

function [bestSig, bestLab, onsets, behDat, fs, figR] = extractSubstrate(od, id)
    assert(isfield(od,'data')&&isfield(od,'labels')&&isfield(od,'fs'), 'missing data/labels/fs');
    fs = od.fs;
    labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
    % bestMac (OB) by label; never by index
    assert(isfield(od,'bestMac')&&~isempty(od.bestMac), 'bestMac field missing');
    bestLab = char(string(od.bestMac));
    bi = find(strcmp(bestLab, labs), 1);
    assert(~isempty(bi), 'bestMac label "%s" not found in labels', bestLab);
    bestSig = double(od.data(bi,:));
    % onsets: cue = one cued sniff/trial (behDat.finalOnset)
    assert(isfield(od,'behDat')&&istable(od.behDat)&&ismember('finalOnset',od.behDat.Properties.VariableNames), ...
        'behDat.finalOnset missing');
    bd = od.behDat; keep = true(height(bd),1);
    if ismember('sniffLabel', bd.Properties.VariableNames)
        sl = string(bd.sniffLabel); if any(sl=="cued"), keep = keep & (sl=="cued"); end
    end
    behDat = bd(keep,:); onsets = double(bd.finalOnset(keep));
    figR = ''; if isfield(od,'figs'), figR = char(string(od.figs)); end
end

function r = summaryRow(id, bestLab, N, retention, pp)
    r = {string(id), string(bestLab), N, retention, nanmean(pp.hasPeak), ...
         nanmed(pp.gammaPeakFrequency), nanmed(pp.gammaPeakTime), ...
         nanmed(pp.peakBurstLength), nanmed(pp.ridgeBurstLength), ...
         nanmed(pp.peakPhaseOffsetWide - pp.peakPhaseOnsetWide), ...
         nanmed(pp.ridgePhaseOffsetWide - pp.ridgePhaseOnsetWide)};
end
function v = summaryVars()
    v = {'session','bestMac','nTrials','retention','fracHasPeak', ...
         'medGammaPeakFreq','medGammaPeakTime','medPeakBurstLen','medRidgeBurstLen', ...
         'medPeakPhaseWideWin','medRidgePhaseWideWin'};
end

function v = getf(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function v = nanmed(x),  x=x(~isnan(x)); if isempty(x), v=NaN; else, v=median(x); end, end
function v = nanmean(x), x=x(~isnan(x)); if isempty(x), v=NaN; else, v=mean(x);   end, end
function v = getenvOr(n,d), v=getenv(n); if isempty(v), v=d; end, v=strtrim(v); end
