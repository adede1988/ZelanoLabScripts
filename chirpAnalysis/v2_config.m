function C = v2_config(overrides)
% V2_CONFIG  Central parameters for chirpAnalysisV2 (spec chirpAnalysis_2.md).
%   C = v2_config()             -> defaults
%   C = v2_config(struct(...))  -> defaults with named fields overridden
%
%   One place for every lever. Signal path is task-agnostic; cueTask is the only current target.

    C = struct();

    % ---- sampling / epoching (spec 2.3) ----
    C.fs        = 500;             % Hz (all finals)
    C.epochWin  = [-1.0 3.0];      % s around behDat.finalOnset (core analysis window)
    C.padSec    = 1.5;             % s padding each side for edge-free FIR/hilbert/FASLT (spec 4)
    C.task      = 'cueTask';

    % ---- noise QC (spec 2.4; cue_noise_trials.m) ----
    C.noise.epWin = [-1.0 3.0];    % window scanned for sharp deflections (= analysis window)
    C.noise.K     = 10;            % robust-z threshold; calibrated (>=80% retained dataset-wide)
    C.noise.winMs = 10;            % sliding max-min window (ms)

    % ---- broadband gamma (phase/envelope substrate; spec 4) ----
    C.bbBand    = [15 75];         % Hz zero-phase FIR (firws); no in-band notch (58 ceiling not needed here)
    C.bbFiltType= 'firws';

    % ---- FASLT TFR (spec 4; reuses chirp_faslt_bank/apply) ----
    C.faslt.range = [20 70];       % Hz (pad beyond 25-58 for context)
    C.faslt.Nf    = 51;            % 1 Hz step over 20-70
    C.faslt.c1    = 3;             % base cycles
    C.faslt.ord   = [3 30];        % fractional adaptive order interval
    C.faslt.mult  = 1;             % multiplicative superresolution

    % ---- baseline z-score (spec 4; myChanZscore) ----
    C.baselineMs  = [-1000 -500];  % baseline period for per-frequency bootstrap z-score

    % ---- ridge extraction (spec 4) ----
    C.ridge.band     = [25 58];    % Hz search band for ridges
    C.ridge.penalty  = 1.0;        % tfridge frequency-jump penalty
    C.ridge.fwhmPeel = true;       % Gaussian-FWHM peel of primary before 2nd ridge
    C.ridge.minPeakZ = 0.5;        % only peel where primary z-power exceeds this (skip noise floor)

    % ---- power/phase continuity (spec 4) ----
    C.pp.smoothWin      = 50;      % samples, gaussian smoothing of ridge power (smoothdata)
    C.pp.peakWinMs      = [0 700]; % search window for the gamma power peak
    C.pp.burstZ         = 2;       % z threshold for burst onset/offset at gammaPeakFrequency
    C.pp.narrowBWhz     = 5;       % narrowband width (+-2.5) around ridge freq
    C.pp.phaseAvgSamp   = 5;       % running average (10 ms) of the phase difference before thresholding
    C.pp.threshNarrow   = pi/4;    % tighter phase-divergence threshold
    C.pp.threshWide     = pi/2;    % looser phase-divergence threshold

    % ---- figures (spec 2.5) ----
    C.figSub    = 'singleTrialSpectrograms';   % per-session subfolder
    C.nFigTrials= 12;              % # example single-trial TF plots to save

    % ---- roots (env overrides for the E: lab session that cannot see R:) ----
    C.dataRoot  = getenvOr('CHIRP_DATAROOT', 'E:\Lab_Common');  % cue finals on E:
    C.figRootE  = getenvOr('CHIRP_FIGROOT_E', 'E:\Lab_Common\Adam\Dupi_processing'); % E: fig mirror
    C.freshCutoff = datetime(2026,6,29,0,0,0);  % use the June-29 cue reruns

    if nargin >= 1 && ~isempty(overrides) && isstruct(overrides)
        C = mergeStruct(C, overrides);
    end
end

function v = getenvOr(name, def)
    v = getenv(name); if isempty(v), v = def; end
    v = strtrim(v);
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
