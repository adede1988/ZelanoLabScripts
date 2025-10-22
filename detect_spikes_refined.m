function [flags, info] = detect_spikes_refined(data, threshold, varargin)
% DETECT_SPIKES_REFINED
% Detect spike-like transients in C x T x N data using moving (max-min),
% then suppress likely gamma oscillations using:
%   (A) isolation/prominence filtering on the moving-range signal
%   (B) optional gamma-envelope persistence mask via bandpass + Hilbert
%
% Usage (simple):
%   flags = detect_spikes_refined(data, 50);
%
% Usage (with gamma mask):
%   flags = detect_spikes_refined(data, 50, 'Fs',1000,'UseGammaMask',true);
%
% Key parameters (name-value):
%   'WinSize'        (odd, default 5)   - window length for moving max-min
%   'RequireFull'    (default true)     - ignore edge samples where full centered window not available
%   'IsolationRadius'(default 5)        - radius (samples) around a detection used to require isolation
%   'MinPromFrac'    (default 0.3)      - center range must exceed max(neighbors) by this fraction
%   'CullRunsLen'    (default 6)        - remove detections that fall in runs >= this many consecutive samples
%   'UseGammaMask'   (default false)    - enable gamma-masking
%   'Fs'             (Hz, required if UseGammaMask=true)
%   'GammaBand'      ([40 90] default)  - band for gamma (Hz)
%   'MinCycles'      (default 3)        - require >= this many cycles of high envelope
%   'EnvThreshK'     (default 2.5)      - envelope threshold = median + k*MAD (per chan/trial)
%   'ZeroCrossMin'   (default 4)        - require >= this many zero-crossings in the window
%
% Outputs:
%   flags  - C x T x N binary (1=spike)
%   info   - struct with intermediates: rangeVals, initialDetections, afterIso, gammaMask

% ---------- Parse inputs ----------
p = inputParser;
p.addParameter('WinSize', 5, @(x)isscalar(x) && mod(x,2)==1);
p.addParameter('RequireFull', true, @(x)islogical(x));
p.addParameter('IsolationRadius', 5, @(x)isscalar(x) && x>=0);
p.addParameter('MinPromFrac', 0.3, @(x)isscalar(x) && x>=0);
p.addParameter('CullRunsLen', 6, @(x)isscalar(x) && x>=1);

p.addParameter('UseGammaMask', false, @(x)islogical(x));
p.addParameter('Fs', [], @(x)isscalar(x) && x>0);
p.addParameter('GammaBand', [40 90], @(x)isnumeric(x) && numel(x)==2 && x(1)>0 && x(2)>x(1));
p.addParameter('MinCycles', 3, @(x)isscalar(x) && x>=1);
p.addParameter('EnvThreshK', 2.5, @(x)isscalar(x) && x>0);
p.addParameter('ZeroCrossMin', 4, @(x)isscalar(x) && x>=0);

p.parse(varargin{:});
S = p.Results;

[C,T,N] = size(data);

% Accept scalar or per-channel thresholds
if isscalar(threshold)
    thr = threshold;
else
    threshold = threshold(:);
    if numel(threshold) ~= C
        error('Vector threshold must be Cx1 (one per channel).');
    end
end

% ---------- Step 1: moving (max-min) ----------
half = (S.WinSize-1)/2;
rangeVals = movmax(data, [half half], 2) - movmin(data, [half half], 2);

% Initial detections (before refinement)
if isscalar(threshold)
    det0 = rangeVals > thr;
else
    thrMat = repmat(reshape(threshold, [C 1 1]), 1, T, N);
    det0   = rangeVals > thrMat;
end
if S.RequireFull && half>0
    det0(:,1:half,:) = false;
    det0(:,T-half+1:T,:) = false;
end

% ---------- Step 2A: Isolation & prominence on range ----------
% Keep only local maxima of range and require "prominence" vs neighbors within radius.
% Also cull long runs of consecutive detections (typical for oscillations).

% Local-max along time (strict)
isLocalMax = false(C,T,N);
if T>=3
    isLocalMax(:,2:T-1,:) = rangeVals(:,2:T-1,:) > rangeVals(:,1:T-2,:) & ...
                            rangeVals(:,2:T-1,:) > rangeVals(:,3:T,:);
end

% Neighborhood max excluding center
rad = S.IsolationRadius;
rangePad = padarray(rangeVals, [0 rad 0], -inf, 'both');  % pad in time
neighMax = zeros(C,T,N);
for offset = [-rad:-1 1:rad]
    neighMax = max(neighMax, rangePad(:, (1+rad+offset):(T+rad+offset), :));
end

promOK = (rangeVals >= (1+S.MinPromFrac).*neighMax);  % center >> neighbors
isoCandidates = det0 & isLocalMax & promOK;

% Cull long runs: compute run-lengths along time of det0 and remove positions in runs >= CullRunsLen
detRuns = det0;
for ch=1:C
  for tr=1:N
    v = squeeze(detRuns(ch,:,tr));
    if any(v)
        % find runs
        dv = diff([false; v(:); false]);
        runStarts = find(dv==1);
        runEnds   = find(dv==-1)-1;
        runLens   = runEnds - runStarts + 1;
        longIdx = find(runLens >= S.CullRunsLen);
        for k=1:numel(longIdx)
            idx = runStarts(longIdx(k)):runEnds(longIdx(k));
            % knock out all detections in long runs (typical oscillations)
            isoCandidates(ch, idx, tr) = false;
        end
    end
  end
end

afterIso = isoCandidates;

% ---------- Step 2B (optional): Gamma-burst mask ----------
% Identify sustained gamma: high Hilbert envelope for >= MinCycles cycles AND multiple zero-crossings.
if S.UseGammaMask
    if isempty(S.Fs)
        error('Fs is required when UseGammaMask=true.');
    end
    Fs = S.Fs;
    % Bandpass
    [b,a] = butter(4, S.GammaBand/(Fs/2), 'bandpass');
    % gammaSig = filtfilt(b,a, data, [], 2);  % C x T x N

     gammaSig = arrayfun(@(x) ...
            filtfilt(b,a, squeeze(data(x,:,:))),...
                1:C, 'uniformoutput', false); 
    gammaSig = cat(3, gammaSig{:}); 
    gammaSig = permute(gammaSig, [3,1,2]); 
    % Envelope threshold per channel/trial: median + k*MAD
    env = arrayfun(@(x) abs(hilbert(squeeze(gammaSig(x,:,:)))),...
        1:C, 'uniformoutput', false);    % C x T x N
    env = cat(3, env{:}); 
    env = permute(env, [3,1,2]); 
    gammaMask = false(C,T,N);

    minDur = ceil(S.MinCycles * Fs / mean(S.GammaBand));  % ~ cycles worth of samples

    for ch=1:C
      for tr=1:N
        e = squeeze(env(ch,:,tr));
        medE = median(e);
        madE = mad(e,1);
        eThr = medE + S.EnvThreshK * madE;

        % candidate high-envelope regions
        hi = e > eThr;
        if ~any(hi), continue; end

        % zero-crossings count in bandpassed signal
        g = squeeze(gammaSig(ch,:,tr));
        zc = sum(abs(diff(sign(g)))==2); %#ok<NASGU> % not used globally; we need per-window below

        % consolidate runs of high envelope and check persistence + zc
        dv = diff([false; hi(:); false]);
        sIdx = find(dv==1);
        eIdx = find(dv==-1)-1;
        for k=1:numel(sIdx)
            seg = sIdx(k):eIdx(k);
            if numel(seg) >= minDur
                % Require enough zero-crossings within this segment (periodicity)
                zcSeg = sum(abs(diff(sign(g(seg))))==2);
                if zcSeg >= S.ZeroCrossMin
                    gammaMask(ch, seg, tr) = true;
                end
            end
        end
      end
    end

    % Remove spike detections that occur inside gamma-mask spans
    afterGamma = afterIso & ~gammaMask;
else
    gammaMask = false(C,T,N);
    afterGamma = afterIso;
end

% ---------- Final flags ----------
flags = double(afterGamma);

% ---------- Info outputs ----------
info = struct();
info.rangeVals          = rangeVals;
info.initialDetections  = det0;
info.afterIsolation     = afterIso;
info.gammaMask          = gammaMask;

end
