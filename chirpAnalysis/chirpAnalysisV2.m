function chirpAnalysisV2(task, sessFilter, opts)
% CHIRPANALYSISV2  Driver for the OB sniff-locked gamma ridge + phase/power continuity analysis,
%   generalized across tasks and split by trial-type category (spec chirpAnalysis_4.md).
%   chirpAnalysisV2()                    % all ready tasks (cueTask, O15) all sessions
%   chirpAnalysisV2('O15')               % one task
%   chirpAnalysisV2('cueTask', {'<id>'}) % one/few sessions (profiling)
%   chirpAnalysisV2('cueTask', [], struct('saveFig',true,'saveMat',false,'maxSessions',1))
%
%   Per session: load final -> bestMac + behDat.finalOnset -> cue_noise_trials QC -> v2_tfr
%   (FASLT + baseline z; O15 long-window guard) -> v2_ridges (primary + peel -> secondary) ->
%   v2_powerphase (gammaPeak, peak/ridge burst, phase continuity, burstTruncated) ->
%   v2_peaklocked (peakLockedRidge + per-trial peak-aligned TFR). Writes ridgeInfo +
%   powerPhaseContinuity into outDat (all trials) and re-saves the final. Then, PER TRIAL-TYPE
%   CATEGORY (cue: hit/cr; O15: start/free/confirm; thresh: air/low/med), computes subject means
%   and figures (filenames tagged <task>_<cat>) and a per-subject aggregate for the group step.
%   ALL helpers called from here (spec 2.6). Groups = outDat.type (Dupi vs OBE).

    if nargin<1 || isempty(task), task = 'all'; end
    if nargin<2, sessFilter = []; end
    if nargin<3, opts = struct(); end
    setup_chirpAnalysis_paths(false);
    here = fileparts(mfilename('fullpath')); repo = fileparts(here);
    addpath(fullfile(repo,'cueAnalysis'));            % cue_noise_trials.m
    C = v2_config();
    saveFig = getf(opts,'saveFig',true); saveMat = getf(opts,'saveMat',true);

    tasks = cellstr(task);
    if strcmpi(task,'all'), tasks = {'cueTask','O15'}; end   % threshTask deferred (no bestMac yet)
    outDir = getenvOr('CHIRP_V2OUT','E:\chirpV2out'); if ~isfolder(outDir), mkdir(outDir); end

    for ti = 1:numel(tasks)
        tk = tasks{ti}; tc = v2_taskconfig(tk);
        files = listFinals(C, tc.suffix);
        if ~isempty(sessFilter), files = files(ismember({files.sessID}, cellstr(string(sessFilter)))); end
        if isfield(opts,'maxSessions') && numel(files) > opts.maxSessions, files = files(1:opts.maxSessions); end
        aggDir = fullfile(outDir, ['agg_' tk]); if ~isfolder(aggDir), mkdir(aggDir); end
        fprintf('\n#### TASK %s : %d session(s) ####\n', tk, numel(files));

        for si = 1:numel(files)
            id = files(si).sessID; fp = files(si).path; t0 = tic;
            fprintf('\n==== [%s] %d/%d  %s ====\n', tk, si, numel(files), id);
            try
                s = load(fp); fn = fieldnames(s); vname = fn{1}; od = s.(vname); clear s;
                [ok,bestSig,bestLab,onsets,catLab,fs,grp,endLimit,msg] = extractSubstrate(od, tc);
                if ~ok, fprintf('  SKIP: %s\n', msg); continue; end
                grp = v2_grouplabel(grp, id);           % Dupi -> DupiS1/S2/S3; OBE stays OBE
                N = numel(onsets);

                nq = cue_noise_trials(bestSig, fs, onsets, C.noise.epWin, C.noise.K, C.noise.winMs);
                good = nq.ok & ~nq.noisy;
                fprintf('  bestMac=%s N=%d retained=%.0f%% grp=%s cats:', bestLab, N, 100*mean(good), grp);

                T   = v2_tfr(bestSig, fs, onsets, good, C, endLimit);
                ridgeInfo = v2_ridges(T, C);
                ridgeInfo.bestMac = bestLab; ridgeInfo.task = tk; ridgeInfo.noiseK = C.noise.K;
                ridgeInfo.noisyMask = nq.noisy(:); ridgeInfo.okMask = nq.ok(:);
                pp  = v2_powerphase(T, ridgeInfo, C);
                [plr, ~, relFreq, relTimeMs, centStack] = v2_peaklocked(T, pp, C);
                ridgeInfo.peakLockedRidge = plr;

                od.ridgeInfo = ridgeInfo; od.powerPhaseContinuity = pp;
                if saveMat, tmp = struct(); tmp.(vname) = od; save(fp,'-struct','tmp','-v7.3'); end

                % peak-centered phase-consistency matrices (once per session)
                [Mpk, phaseRelT] = v2_peakcenter(pp.peakPhaseConsistency,  pp.tMs, pp.gammaPeakTime, 1000, fs);
                [Mrg, ~        ] = v2_peakcenter(pp.ridgePhaseConsistency, pp.tMs, pp.gammaPeakTime, 1000, fs);

                figDirE = fullfile(C.figRootE, id, C.figSub);
                agg = struct('sessID',id,'group',grp,'task',tk,'bestMac',bestLab, ...
                    'freqs',T.freqs,'tMs',T.tMs,'relFreq',relFreq,'relTimeMs',relTimeMs,'phaseRelT',phaseRelT);
                byCat = struct('cat',{},'meanTFR',{},'meanPeakTFR',{},'peakConsMean',{},'ridgeConsMean',{}, ...
                    'peakBurst',{},'ridgeBurst',{},'medPeakBurst',{},'medRidgeBurst',{},'nTrials',{});
                catStr = "";
                for ci = 1:numel(tc.cats)
                    cat = tc.cats{ci};
                    gc = find(good & (catLab==string(cat)));
                    if isempty(gc), continue; end
                    catStr = catStr + sprintf(' %s=%d', cat, numel(gc));
                    tag = sprintf('%s_%s', tk, cat); tagT = strrep(tag,'_','\_');
                    meanTFR     = mean(T.zTFR(:,:,gc), 3, 'omitnan');
                    meanPeakTFR = mean(centStack(:,:,gc), 3, 'omitnan');
                    if saveFig
                        % single-trial examples (spread across the category's trials)
                        pick = gc(unique(round(linspace(1, numel(gc), min(6,numel(gc))))));
                        for j = pick(:)'
                            v2_fig_singletrial(T.zTFR(:,:,j), T.freqs, T.tMs, ...
                                ridgeInfo.primaryRidge.f(j,:), ridgeInfo.secondaryRidge.f(j,:), ...
                                fullfile(figDirE, sprintf('sub-%s_%s_ch-%s_trial-%03d.png', id, tag, bestLab, j)), ...
                                sprintf('%s  %s  trial %d', strrep(id,'_','\_'), tagT, j), C);
                        end
                        v2_fig_meantfr(meanTFR, T.freqs, T.tMs, ...
                            fullfile(figDirE, sprintf('meanTFR_%s_%s.png', tag, id)), ...
                            sprintf('%s  %s  mean TFR (n=%d)', strrep(id,'_','\_'), tagT, numel(gc)), C);
                        v2_fig_phasecons(Mpk, Mrg, gc, phaseRelT, ...
                            fullfile(figDirE, sprintf('phaseConsistency_%s_%s.png', tag, id)), ...
                            sprintf('%s  %s  phase continuity (n=%d)', strrep(id,'_','\_'), tagT, numel(gc)), C);
                        v2_fig_bursthist(pp.peakBurstLength(gc), pp.ridgeBurstLength(gc), ...
                            fullfile(figDirE, sprintf('burstLength_%s_%s.png', tag, id)), ...
                            sprintf('%s  %s  burst length', strrep(id,'_','\_'), tagT));
                        v2_fig_peaktfr(meanPeakTFR, relFreq, relTimeMs, ...
                            fullfile(figDirE, sprintf('peakAlignedTFR_%s_%s.png', tag, id)), ...
                            sprintf('%s  %s  peak-aligned mean TFR', strrep(id,'_','\_'), tagT));
                    end
                    byCat(end+1) = struct('cat',cat,'meanTFR',meanTFR,'meanPeakTFR',meanPeakTFR, ...
                        'peakConsMean',mean(abs(Mpk(gc,:)),1,'omitnan'),'ridgeConsMean',mean(abs(Mrg(gc,:)),1,'omitnan'), ...
                        'peakBurst',pp.peakBurstLength(gc(~isnan(pp.peakBurstLength(gc)))), ...
                        'ridgeBurst',pp.ridgeBurstLength(gc(~isnan(pp.ridgeBurstLength(gc)))), ...
                        'medPeakBurst',nanmed(pp.peakBurstLength(gc)),'medRidgeBurst',nanmed(pp.ridgeBurstLength(gc)), ...
                        'nTrials',numel(gc)); %#ok<AGROW>
                end
                agg.byCat = byCat;
                save(fullfile(aggDir,[id '_agg.mat']),'agg','-v7.3');
                fprintf('%s | done %.1fs\n', catStr, toc(t0));
                clear T ridgeInfo pp od centStack Mpk Mrg;
            catch ME
                fprintf('  *** SESSION FAILED (%s): %s\n', id, ME.message);
                fprintf('      %s\n', getReport(ME,'basic','hyperlinks','off'));
            end
        end

        if getf(opts,'doGroup', numel(files)>1)
            try
                v2_group(aggDir, fullfile(C.figRootE,'groupStatFigs'), outDir, C, tk, tc.cats);
            catch GE
                fprintf('group step failed (%s): %s\n', tk, GE.message);
            end
        end
    end
    fprintf('\nchirpAnalysisV2 complete.\n');
end

% =================== session listing / substrate ===================
function files = listFinals(C, suffix)
    d = dir(fullfile(C.dataRoot,'**',['*_' suffix '.mat']));   % case-insensitive on Windows
    files = struct('sessID',{},'path',{});
    for k = 1:numel(d)
        if datetime(d(k).datenum,'ConvertFrom','datenum') < C.freshCutoff, continue; end
        [~,nm] = fileparts(d(k).name);
        id = regexprep(nm, ['_' suffix '$'], '', 'ignorecase');
        files(end+1) = struct('sessID',id,'path',fullfile(d(k).folder,d(k).name)); %#ok<AGROW>
    end
end

function [ok,bestSig,bestLab,onsets,catLab,fs,grp,endLimit,msg] = extractSubstrate(od, tc)
    ok=false; bestSig=[]; bestLab=''; onsets=[]; catLab=strings(0); fs=[]; grp='Dupi'; endLimit=[]; msg='';
    if ~(isfield(od,'data')&&isfield(od,'labels')&&isfield(od,'fs')), msg='missing data/labels/fs'; return; end
    fs = od.fs;
    labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
    if ~(isfield(od,'bestMac')&&~isempty(od.bestMac)), msg='no bestMac (skip)'; return; end
    bestLab = char(string(od.bestMac));
    bi = find(strcmp(bestLab, labs), 1);
    if isempty(bi), msg=sprintf('bestMac "%s" not in labels', bestLab); return; end
    bestSig = double(od.data(bi,:));
    if ~(isfield(od,'behDat')&&istable(od.behDat)&&ismember('finalOnset',od.behDat.Properties.VariableNames)), msg='no behDat.finalOnset'; return; end
    bd = od.behDat;
    if ~ismember(tc.catCol, bd.Properties.VariableNames), msg=sprintf('behDat missing %s', tc.catCol); return; end
    onsets = double(bd.finalOnset);          % keep ALL rows (aligned to behDat); invalid -> NaN downstream
    catLab = string(bd.(tc.catCol));
    if isfield(od,'type')&&~isempty(od.type), grp = char(string(od.type)); end
    if tc.guard, endLimit = nextOnsetLimit(onsets); end
    ok = true;
end

function el = nextOnsetLimit(onsets)
    el = inf(numel(onsets),1);
    for i = 1:numel(onsets)
        nxt = onsets(onsets > onsets(i));
        if ~isempty(nxt), el(i) = min(nxt); end
    end
end

% =================== small helpers ===================
function v = getf(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function v = nanmed(x),  x=x(~isnan(x)); if isempty(x), v=NaN; else, v=median(x); end, end
function v = getenvOr(n,d), v=getenv(n); if isempty(v), v=d; end, v=strtrim(v); end
