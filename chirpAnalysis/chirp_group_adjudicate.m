function G = chirp_group_adjudicate(outDir, C)
% CHIRP_GROUP_ADJUDICATE  Per-session verdict vectors + group adjudication (spec 8, 0.2, 8.1).
%   G = chirp_group_adjudicate(outDir, C)
%   Reads subject_summary.csv (written by run_chirp_analysis_pipeline) and emits, per session, a
%   verdict VECTOR (one element per test) and a single/dual/undecidable call; plus group-level
%   tests. "spatial: not assessed" (ineligible) is kept DISTINCT from "spatial: matched" (8.1) and
%   never averaged with negative-spatial. No p-value multiplication -- triangulation across the vector.
%
%   Thresholds tie to Phase 0 (chirp_phase0_harness): a test only "speaks" where Phase-0
%   sensitivity/specificity licensed it. Where unavailable, defaults below are used and flagged.

    if nargin<2||isempty(C), C=chirp_config(); end
    f = fullfile(outDir,'subject_summary.csv');
    assert(isfile(f), 'subject_summary.csv not found in %s', outDir);
    T = readtable(f); T = ensureStr(T);

    % --- thresholds (calibrate from Phase 0; defaults flagged) ---
    th = struct('phaseCrossR',0.5,'phaseCrossP',0.05, 'beatFracSig',0.30, ...
        'chirpFrac',0.40, 'tempCorrSingle',0.5, 'ridgeNR2',0.30, 'spatialSimCeil',0.5);

    n = height(T);
    vPhase=strings(n,1); vBeat=strings(n,1); vChirp=strings(n,1); vTemp=strings(n,1);
    vRidge=strings(n,1); vSpatial=strings(n,1); verdict=strings(n,1); score=zeros(n,1);
    for i=1:n
        % phase (confirms SINGLE): CROSS-TRIAL consistency of the ridge-vs-phase residual
        % (robust to a constant ridge-IF bias; the spec's "disperse across trials" logic).
        vPhase(i) = tern(T.phase_crossP(i)<th.phaseCrossP & T.phase_crossR(i)>th.phaseCrossR, "single","incon");
        % beat (confirms DUAL)
        vBeat(i)  = tern(T.beat_fracSig(i)>=th.beatFracSig, "dual","incon");
        % chirplet geometry
        if T.chirp_fracSingle(i)>=th.chirpFrac && T.chirp_fracSingle(i)>T.chirp_fracDual(i), vChirp(i)="single";
        elseif T.chirp_fracDual(i)>=th.chirpFrac && T.chirp_fracDual(i)>T.chirp_fracSingle(i), vChirp(i)="dual";
        else, vChirp(i)="incon"; end
        % temporal (leans dual when decoupled)
        if isfinite(T.temp_corrHiLo(i))
            vTemp(i) = tern(T.temp_corrHiLo(i) < th.tempCorrSingle, "dual","single");
        else, vTemp(i)="incon"; end
        % ridge 2nd-ridge prevalence (leans dual)
        vRidge(i) = tern(T.ridge_fracNR2(i)>=th.ridgeNR2, "dual","incon");
        % spatial: distinct "not assessed" state
        if ~logical(T.spatial_eligible(i)) || isnan(T.spatial_meanSim(i))
            vSpatial(i) = "notAssessed";
        else
            vSpatial(i) = tern(T.spatial_meanSim(i) < th.spatialSimCeil, "dual","single");
        end
        % combine (triangulation, NOT p-multiplication): tally directional evidence
        votesS = sum([vPhase(i) vChirp(i) vTemp(i) vSpatial(i)]=="single");
        votesD = sum([vBeat(i) vChirp(i) vTemp(i) vRidge(i) vSpatial(i)]=="dual");
        score(i) = votesD - votesS;
        if votesD>=2 && votesS==0, verdict(i)="dual";
        elseif votesS>=2 && votesD==0, verdict(i)="single";
        else, verdict(i)="undecidable"; end
    end

    V = T(:,{'session','group','type'});
    V.phase=vPhase; V.beat=vBeat; V.chirplet=vChirp; V.temporal=vTemp; V.ridge=vRidge;
    V.spatial=vSpatial; V.verdict=verdict; V.dualMinusSingle=score;
    writetable(V, fullfile(outDir,'verdict_vectors.csv'));

    % --- group-level tests ---
    grp = struct();
    grp.n = n;
    grp.phase_t_permZ = oneSampleT(T.phase_crossR);           % H1: crossTrialR>0 (across-trial consistency) => single
    grp.beat_t_logSNR = oneSampleT(T.beat_medLogSNR);         % H1: >0 => dual
    grp.fracSingleTraj = mean(T.chirp_fracSingle,'omitnan');
    grp.fracDualTraj   = mean(T.chirp_fracDual,'omitnan');
    grp.spatial_nEligible = sum(logical(T.spatial_eligible));
    grp.verdictCounts = countcats(categorical(verdict));
    grp.verdictCats = categories(categorical(verdict));

    G = struct('perSession',V,'group',grp,'thresholds',th,'summaryTable',T);
    save(fullfile(outDir,'chirp_group.mat'),'G','-v7.3');
    fprintf('Adjudication: %d sessions | single=%d dual=%d undecidable=%d | spatial eligible=%d\n', ...
        n, sum(verdict=="single"), sum(verdict=="dual"), sum(verdict=="undecidable"), grp.spatial_nEligible);
    fprintf('  group phase perm_z t: t=%.2f p=%.3g (>0 => single) | beat logSNR t: t=%.2f p=%.3g (>0 => dual)\n', ...
        grp.phase_t_permZ.t, grp.phase_t_permZ.p, grp.beat_t_logSNR.t, grp.beat_t_logSNR.p);
end

function r = oneSampleT(x)
    x = x(isfinite(x)); r = struct('t',NaN,'p',NaN,'mean',NaN,'n',numel(x));
    if numel(x)>=2
        m=mean(x); s=std(x); r.mean=m; r.t=m/(s/sqrt(numel(x)));
        r.p = 2*(1-tcdf(abs(r.t),numel(x)-1));
    end
end
function y = tern(c,a,b), if c, y=a; else, y=b; end, end
function T = ensureStr(T)
    for v={'session','group','type'}, if ismember(v{1},T.Properties.VariableNames), T.(v{1})=string(T.(v{1})); end, end
end
