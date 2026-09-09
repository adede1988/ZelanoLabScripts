function C = chirp_config(overrides)
% CHIRP_CONFIG  Central default parameters for the chirp (one-vs-two oscillator) analysis.
%   C = chirp_config()            -> defaults (spec Appendix + design-panel/critic decisions)
%   C = chirp_config(struct(...)) -> defaults with the named fields overridden
%
%   Every analysis function takes this struct so parameters live in ONE place (tuned in
%   Phase 0). Output/data roots honor environment overrides so the SSH (no-R:) lab session
%   writes to E:.  See chirpAnalysisStatus.md sec "DESIGN-PANEL FIXES" for why each lever is
%   set the way it is. Tune in Phase 0; do not scatter magic numbers in the test functions.

    C = struct();

    % ---- sampling / epoching ----
    C.fs          = 500;            % Hz (fixed by the pipeline)
    C.epochWin    = [-1.0 3.0];     % s around behDat.finalOnset (analysis window; D3)
    C.padSec      = 1.0;            % s extra each side for edge-free FIR/TFR (>= max FASLT wavelet
                                    %   half-support in the 25-58 band ~0.73 s; chirp_tfr_faslt asserts).
    C.task        = 'cueTask';

    % ---- broadband gamma (phase & envelope substrate; ALWAYS hilbert of this, NOT FASLT) ----
    C.bbBand      = [25 58];        % Hz, zero-phase FIR (firws); 58 ceiling clears 60 Hz line
    C.bbFiltType  = 'firws';

    % ---- FASLT time-frequency (TFR substrate: spectrogram display + ridge ONLY) ----
    C.faslt.range = [20 70];        % Hz (pad beyond 25-58 for context)
    C.faslt.Nf    = 51;             % => 1 Hz step over 20-70
    C.faslt.c1    = 3;              % base cycles
    C.faslt.ord   = [3 30];         % fractional adaptive order interval (linear over freq; 25 Hz~5.7, 58 Hz~23)
    C.faslt.mult  = 1;             % multiplicative superresolution
    %   NOTE (Phase-0 lever): high ohi -> long wavelet -> poorer time res for a fast sweep. If the
    %   ridge cannot track the synthetic nonlinear sweep, lower ohi (e.g. [2 10]).

    % ---- ridge (tfridge on FASLT power, 25-58 band) ----
    C.ridge.band       = [25 58];   % Hz search band for ridges
    C.ridge.penalty    = 1.0;       % tfridge frequency-jump penalty (CALIBRATE in Phase 0)
    C.ridge.maxRidges  = 2;         % primary + 1 peeled
    C.ridge.nSurr      = 60;        % surrogates per sampled trial for the per-CHANNEL ridge-2 null (F6)
    C.ridge.nNullTrials= 3;         % # sampled trials pooled for the per-channel ridge-2 null ratio
    C.ridge.peelHalfHz = 4;         % +-Hz removed around the primary ridge (full notch ~8 Hz < min Df)
    C.ridge.anchorClamp= [27 56];   % Hz clamp for f_hi/f_lo anchors (subject-avg TFR; ANNOTATION ONLY)
    C.ridge.transHiPct = 20;        % anchor f_hi read at ~20th pctile transition time
    C.ridge.transLoPct = 80;        % anchor f_lo read at ~80th pctile transition time
    C.ridge.fhatSmoothMs = 75;      % movmedian smoothing of fhat before transition detection
    C.ridge.transEdgeHz  = 1;       % transition = fhat from (f_hi-1) down to within 1 Hz of f_lo plateau
    %   CRITIC FIX (high): per-trial transitionWin is derived from THIS trial's own smoothed fhat
    %   extrema, NEVER from the inward-biased subject-average anchors (avoids a selection bias).

    % ---- phase continuity (PRIMARY single-channel test) ----
    C.phase.subWinMs    = 150;      % consecutive sub-window length (<< 1-cycle IF-bias accumulation)
    C.phase.minTransSec = 0.30;     % require >=2 sub-windows of transition before a trial contributes
    C.phase.ifSource    = 'ridgePower';  % f-hat = POWER ridge, NEVER dphi/dt (anti-bias rule a) -- asserted
    C.phase.maxRidgeGap = 0.25;     % drop trial if ridge NaN over >25% of its transition
    C.phase.slipThreshRad = pi/2;   % |sub-window resid| past this = a phase slip (localizes beat nulls)
    C.phase.requireResponse = true; % existence-only gate (rule b); never shape-based
    C.phase.minRespZ    = 1.0;      % median in-transition gamma power must exceed 1 SD over baseline
    C.phase.nPerm       = 1000;     % per-session rotation permutation null -> perm_z

    % ---- envelope beating (confirms dual) ----
    C.beat.envDomain   = 'amplitude'; % beat env = abs(hilbert(broadband)); decoup/spatial use power
    C.beat.envPolyOrder= 3;         % low-order poly detrend (remove onset/offset shape)
    C.beat.hpCutoffHz  = [];        % NO highpass by default; if set must be < min Df (<=3), NEVER 10
    C.beat.minCycles   = 2;         % require >=2 beat cycles in the COEXISTENCE window
    C.beat.flankHz     = 3;         % half-width of the flanking comparison band around fBeat

    % ---- chirplet (MPACT MPEM) -- judge GEOMETRY not count (rule c) ----
    C.chirplet.Q        = 3;        % up to 3 atoms so a curved SINGLE sweep can use 2-3 LINKED atoms
    C.chirplet.M        = 64;       % Newton-Raphson refinement resolution (MPACT default; Phase-0 lever)
    C.chirplet.D        = 5;        % decomposition depth
    C.chirplet.i0       = 1;        % first scale to rotate
    C.chirplet.radix    = 2;
    C.chirplet.mnits    = 5;        % max EM refinement iterations
    C.chirplet.level    = 2;        % MLE difficulty (unused on expectmax branch)
    C.chirplet.useAnalytic = true;  % feed hilbert(seg) (one-sided) to avoid mirror-freq atoms
    C.chirplet.atomKeepBand = [26 57]; % keep atoms with fc in this band (inside FIR band; not [22 62])
    C.chirplet.cTolHzs  = 10;       % |chirp rate| (Hz/s) below which an atom counts as 'flat' (Model-B)
    C.chirplet.crMin    = 2;        % a lone atom must descend faster than this (Hz/s) to count as a sweep
                                    %   (lowered from 5: realistic slow sweeps ~5-10 Hz over 1-2 s = 3-10 Hz/s)
    C.chirplet.dualDf   = 3;        % two coexisting flat atoms must differ in fc by > this (Hz) -> dual
    C.chirplet.dfMin    = 1;        % chained atoms must step DOWN in fc by > this (Hz)  [legacy; see dfTol]
    C.chirplet.dfTol    = 1.5;      % chain edge if consecutive fc is NON-INCREASING within +this Hz
                                    %   (allows a flat plateau / fit jitter -> nonlinear single = one chain)
    C.chirplet.chainGapK= 1.5;      % two atoms chain if tc gap <= this * median(sigma_s)
    C.chirplet.trajThiSingle = 0.66;% trajConnectivity >= this (no parallel pair) -> single
    C.chirplet.minTransSamp  = 75;  % 0.15 s; below this the transition can't host a fit -> included=false
    C.chirplet.NeffMode = 'cycles'; % BIC effective N = #gamma cycles in window, NOT raw samples
    C.chirplet.dBICcorroborateOnly = true; % dBIC never overrides a geometry 'single' call
    C.chirplet.sigMode = 'fixed';   % 'fixed' (fast batch) | 'surrogate' (per-trial; Phase-0 calibration)
    C.chirplet.atomEnergyThr = 0.30;% fixed-mode: keep atoms capturing >30% of signal energy (CALIBRATE
                                    %   Phase 0). 0.3 keeps dominant components, drops residual junk that
                                    %   would spuriously break a single sweep's chain (0.08 was too low).
    C.chirplet.nSurr   = 40;        % surrogate-mode atom-significance surrogates (when sigMode='surrogate')
    C.chirplet.maxWinMs = 600;      % cap fit window to this many ms (bounds MPEM cost; empty = no cap)

    % ---- surrogates (SHARED generator; ridge-2, beat modSNR, chirplet atom) ----
    C.surr.n      = 200;            % per-trial FFT phase-randomized surrogates (batch); 1000 for profiled
    C.surr.pctile = 95;            % significance threshold percentile

    % ---- spatial (across macBP*; per user D2: NO geometry, eligible if >=2 good contacts) ----
    C.spatial.minGoodContacts = 2;  % eligible if >=2 good macBP contacts (>=3 preferred)
    C.spatial.bandHalfHz      = 2.5;% +-Hz around f_hi / f_lo for early-high / late-low band power
    C.spatial.useGeometry     = false; % user waived layout/spacing (D2); contactSpacing = NaN

    % ---- temporal decoupling ----
    C.temporal.bandHalfHz = 2.5;    % +-Hz around f_hi / f_lo for latency envelopes (FIR-hilbert power)
    C.temporal.smoothMs   = 40;     % smooth band-limited power envelope before peak-latency

    % ---- burst-length measurement (D4; user-specified) ----
    C.burst.smoothMs   = 100;       % smoothing of the gamma-power time series before thresholding
    C.burst.peakPct    = 75;        % threshold = this pctile of the CROSS-TRIAL peak distribution
    C.burst.searchMs   = [0 3000];  % search start/stop only within (0, 3000] ms post-onset

    % ---- existence check (optional; NOT shape-based) ----
    C.existence.enforce = false;    % default OFF (a response-free trial just adds phase noise)

    % ---- reproducibility ----
    C.rngSeed = 20260629;           % deterministic per-(session,trial) seeding for surrogates/perm/synth

    % ---- fresh-final filter (D5): use only finals modified after this instant ----
    C.freshCutoff = datetime(2026,6,29,0,0,0);   % all June-29 cue reruns qualify

    % ---- roots (env overrides for the lab SSH session that cannot see R:) ----
    C.dataRoot = getenvOr('CHIRP_DATAROOT', 'E:\Lab_Common');     % where cue finals live on E:
    C.figRoot  = getenvOr('CHIRP_FIGROOT',  '');                  % per-session fig parent; '' -> labPaths.figPath
    C.outDir   = getenvOr('CHIRP_OUTDIR',   '');                  % aggregate CSV dir; '' -> <groupDir>\chirpAnalysis
    C.figSub   = 'singleTrialTF';                                 % per-session fig subfolder (user-requested name)

    % apply overrides (recursive struct merge)
    if nargin >= 1 && ~isempty(overrides) && isstruct(overrides)
        C = mergeStruct(C, overrides);
    end
end

function v = getenvOr(name, def)
    v = getenv(name); if isempty(v), v = def; end
    v = strtrim(v);   % guard against trailing space from cmd `set VAR=val & ...`
end

function A = mergeStruct(A, B)
    fn = fieldnames(B);
    for i = 1:numel(fn)
        if isstruct(B.(fn{i})) && isfield(A, fn{i}) && isstruct(A.(fn{i}))
            A.(fn{i}) = mergeStruct(A.(fn{i}), B.(fn{i}));
        else
            A.(fn{i}) = B.(fn{i});
        end
    end
end
