function rsa_panel()
% Vagal-tone panel screen, spontaneous blocks only, per-pair common support.
% Gate 0 feasibility -> Gate 1 kappa reliability (outcome-blind) -> Gate 2 glimmer.

maxNumCompThreads(3); rng(7);
P      = ohrv_config();
sp     = P.code;
outDir = P.work;

SPONT   = ["audio","focus","naturalFocus"];
DROPSEC = 30;      % adaptation transient dropped from each block
NBOOT   = 300;     % block-bootstrap draws
SEGSEC  = 45;      % bootstrap segment length

%% ---------- pass 1: load, restrict to spontaneous, collect ----------
files = dir(fullfile(outDir, '*_slim.mat'));
Sx = struct('sess', {}, 'subj', {}, 'num', {}, 'T', {}, 'runs', {}, 'fs', {}, 'rrint', {});
allVol = [];
for ii = 1:numel(files)
    L = load(fullfile(files(ii).folder, files(ii).name));
    T = L.T; fs = L.fs; nSamp = L.nSamp; sessID = L.sessID; rrint = L.rrint; clear L
    bf = fullfile(outDir, [sessID '_beats.npz']);
    if ~exist(bf, 'file'), continue; end
    hb = readNPZ(bf);

    keep = T.goodBreath == 1 & isfinite(T.RR_max_min) & T.RR_max_min > 0 ...
         & T.len >= 1.5 & T.len <= 15 & isfinite(T.inhDur) & T.inhDur > 0 & T.inhDur < T.len ...
         & isfinite(T.inhVol) & T.inhVol > 0 ...
         & isfinite(T.finalOnset) & T.finalOnset >= 1 ...
         & (T.finalOnset + round(T.len*fs)) <= nSamp & ~strcmpi(T.noseMouth, "mouth") ...
         & ismember(T.task, SPONT);
    T = T(keep, :);
    if height(T) < 20, continue; end

    % contiguous runs = blocks, first DROPSEC seconds removed
    runs = struct('task', {}, 'nn', {}, 'dur', {}, 'rej', {});
    ub = unique(T.condition);
    for k = 1:numel(ub)
        m  = T.condition == ub(k);
        t0 = min(T.finalOnset(m)) + DROPSEC*fs;
        t1 = max(T.finalOnset(m) + round(T.len(m)*fs));
        bb = hb(hb >= t0 & hb <= t1);
        if numel(bb) < 40, continue; end
        [nn, rej] = cleanNN(diff(bb)/fs);
        if sum(isfinite(nn)) < 30, continue; end
        runs(end+1) = struct('task', T.task(find(m,1)), 'nn', nn, ...
                             'dur', (t1-t0)/fs, 'rej', rej); %#ok<AGROW>
    end
    if isempty(runs), continue; end

    % ---- Gate 0: feasibility, outcome-blind ----
    % Gate on what SURVIVES, not on how much was removed. The rejection rate is
    % a property of the filter's aggressiveness as much as of the data: adding
    % the stage-3 ectopic rule roughly doubled it in clean sessions, which would
    % spuriously fail them under a threshold calibrated to the weaker filter.
    % What actually matters is whether enough clean beats remain and whether the
    % survivors are still outlier-driven (rmssdRatio, the diagnostic that caught
    % the DB artifact).
    totMin = sum([runs.dur]) / 60;
    rejF   = mean([runs.rej]);
    nnAllQ = vertcat(runs.nn); nClean = sum(isfinite(nnAllQ));
    ddQ = []; for q = 1:numel(runs), dq = diff(runs(q).nn); ddQ = [ddQ; dq(isfinite(dq))]; end %#ok<AGROW>
    ratioQ = sqrt(mean(ddQ.^2)) / max(median(abs(ddQ)), eps);
    if height(T) < 60 || totMin < 3 || nClean < 200 || rejF > 0.35 || ratioQ > 3
        fprintf('GATE0 FAIL  %-28s nBr=%3d  min=%4.1f  nClean=%4d  rej=%4.1f%%  ratio=%4.1f\n', ...
            sessID, height(T), totMin, nClean, 100*rejF, ratioQ);
        continue
    end

    tok = split(string(sessID), '_'); sj = tok(4); if sj == "TPB", sj = "TB"; end
    Sx(end+1) = struct('sess', string(sessID), 'subj', sj, 'num', double(tok(5)), ...
                       'T', T, 'runs', runs, 'fs', fs, 'rrint', rrint); %#ok<AGROW>
    allVol = [allVol; T.inhVol]; %#ok<AGROW>
end
VREF = median(allVol); LREF = 4;
fprintf('pass 1: %d sessions with spontaneous data;  vol ref = %.1f\n\n', numel(Sx), VREF);

fprintf('%-30s %-6s %5s %6s %7s %7s  tasks\n', 'session','subj','nBr','nNN','min','rej%');
for i = 1:numel(Sx)
    nn = vertcat(Sx(i).runs.nn); tot = sum([Sx(i).runs.dur])/60;
    fprintf('%-30s %-6s %5d %6d %7.1f %7.1f  %s\n', Sx(i).sess, Sx(i).subj, ...
        height(Sx(i).T), sum(isfinite(nn)), tot, 100*mean([Sx(i).runs.rej]), ...
        strjoin(unique([Sx(i).runs.task]), ','));
end

%% ---------- pass 2: pairs ----------
METS = {'evokedRSA','depthSlope','logRMSSD','logMASD','logMedPV','adjLogRSA','logHF','pNN20', ...
        'logMeanNN','breathRate','rmssdRatio'};
subs = unique([Sx.subj]);
P = table();
for s = 1:numel(subs)
    idx = find([Sx.subj] == subs(s));
    [~, o] = sort([Sx(idx).num]); idx = idx(o);
    for j = 1:numel(idx)-1
        a = Sx(idx(j)); b = Sx(idx(j+1));
        ta = unique([a.runs.task]); tb = unique([b.runs.task]);
        common = intersect(ta, tb);
        if isempty(common), continue; end
        ra = a.runs(ismember([a.runs.task], common));
        rb = b.runs(ismember([b.runs.task], common));
        % duration match: truncate the longer to the shorter total beat count
        na = sum(arrayfun(@(r) sum(isfinite(r.nn)), ra));
        nb = sum(arrayfun(@(r) sum(isfinite(r.nn)), rb));
        cap = min(na, nb);
        ra = truncRuns(ra, cap); rb = truncRuns(rb, cap);

        ma = metrics(ra, a.T, common, VREF, LREF, a.rrint, a.fs);
        mb = metrics(rb, b.T, common, VREF, LREF, b.rrint, b.fs);
        sa = bootSE(ra, a.T, common, VREF, LREF, NBOOT, SEGSEC, METS, a.rrint, a.fs);
        sb = bootSE(rb, b.T, common, VREF, LREF, NBOOT, SEGSEC, METS, b.rrint, b.fs);

        row = table(subs(s), a.num, b.num, cap, strjoin(common, ','), ...
            'VariableNames', {'subj','from','to','nNN','blocks'});
        for k = 1:numel(METS)
            row.(['d_' METS{k}])  = mb.(METS{k}) - ma.(METS{k});
            row.(['se_' METS{k}]) = sqrt(sa.(METS{k})^2 + sb.(METS{k})^2);
        end
        P = [P; row]; %#ok<AGROW>
    end
end
fprintf('\n--- pairs ---\n'); disp(P(:, [1:5, 6:2:6+2*numel(METS)-1]));
writetable(P, fullfile(outDir, 'panel_pairs.csv'));

%% ---------- Gate 1: kappa (outcome-blind) ----------
Q = P(P.from == 1 & P.to == 2, :);       % primary contrast: matched ~1 month
fprintf('\nPRIMARY CONTRAST: session 1 -> 2,  n = %d subjects (%s)\n', ...
    height(Q), strjoin(cellstr(Q.subj)', ', '));
fprintf('\n--- Gate 1: kappa = mean SE^2 / var(delta), outcome-blind ---\n');
keepM = {};
for k = 1:numel(METS)
    d = Q.(['d_' METS{k}]); se = Q.(['se_' METS{k}]);
    kap = mean(se.^2, 'omitnan') / var(d, 'omitnan');
    ok = kap <= 0.5 && isfinite(kap);
    fprintf('  %-12s  var=%8.4f  meanSE=%7.4f  kappa=%6.2f  %s\n', ...
        METS{k}, var(d,'omitnan'), mean(se,'omitnan'), kap, ternS(ok,'RETAIN','drop'));
    % logMeanNN / breathRate are confound companions, rmssdRatio is a QC readout;
    % none of them are candidate vagal metrics, so they never enter the family.
    if ok && ~ismember(METS{k}, {'logMeanNN','breathRate','rmssdRatio'})
        keepM{end+1} = METS{k}; %#ok<AGROW>
    end
end
save(fullfile(outDir, 'panel_state.mat'), 'P', 'Q', 'METS', 'keepM');
fprintf('\nretained for the max statistic: %s\n', ternS(~isempty(keepM), strjoin(keepM, ', '), 'NONE'));
fprintf('PANEL DONE\n');
end

% ================= helpers =================
function [nn, rej] = cleanNN(nn)
% Three-stage artifact rejection.
%
% Stage 2 uses an 11-wide median (spanning a whole breath cycle, so genuine RSA
% averages out of the reference) with a 35% tolerance. A narrower window or
% tighter tolerance rejects real RSA - at slow deep breathing the within-breath
% RR swing routinely exceeds 20%.
%
% Stage 3 is the one that matters and was missing before. An ectopic beat
% produces a PREMATURE short interval followed by a COMPENSATORY long one.
% Each interval individually sits within 35% of the local median, so stages 1-2
% pass it, but their DIFFERENCE is impossible. Because RMSSD squares successive
% differences, a handful of these dominate an entire session: in the previous
% run DB session 1 had RMSSD/MASD = 5.99 (others 1.5-2.1) and a 512 ms maximum
% successive difference, which manufactured a spurious result.
% The threshold is set from a robust scale of the difference series itself, so
% it adapts to each subject's own RSA depth instead of assuming one.
nn = nn(:);
ok = nn > 0.3 & nn < 2.0;                       % stage 1: physiological range
m  = movmedian(nn, 11, 'omitnan');
ok = ok & abs(nn - m) <= 0.35 * m;              % stage 2: local level
nn(~ok) = NaN;

d = diff(nn);
q = median(abs(d), 'omitnan');                  % MASD = robust scale of dRR
if isfinite(q) && q > 0
    thr = 5 * 1.4826 * q;                       % ~5 robust SD
    bad = find(abs(d) > thr);
    for i = bad(:)'                             % drop both intervals spanning the jump
        ok(i) = false; ok(min(i+1, numel(ok))) = false;
    end
end
rej = 1 - mean(ok);
nn(~ok) = NaN;                 % keep positions so diffs across gaps are dropped
end

function r = truncRuns(r, cap)
tot = 0;
for i = 1:numel(r)
    n = sum(isfinite(r(i).nn));
    if tot + n <= cap, tot = tot + n; continue; end
    room = cap - tot; f = find(isfinite(r(i).nn));
    if room > 0, r(i).nn = r(i).nn(1:f(room)); else, r(i).nn = []; end
    r = r(1:i); return
end
end

function M = metrics(runs, T, common, VREF, LREF, rrint, fs)
nn = []; dd = []; hfw = []; wts = [];
for i = 1:numel(runs)
    v = runs(i).nn;
    nn = [nn; v(isfinite(v))]; %#ok<AGROW>
    d  = diff(v); dd = [dd; d(isfinite(d))]; %#ok<AGROW>   % NaN kills cross-gap pairs
    f = v(isfinite(v));
    if numel(f) >= 40
        t = cumsum(f); tu = t(1):0.25:t(end);
        if numel(tu) >= 128
            x = detrend(interp1(t, f, tu, 'pchip'));
            w = min(256, 2^floor(log2(numel(x))));
            [px, fq] = pwelch(x, hamming(w), round(w/2), [], 4);
            hfw(end+1) = trapz(fq, px); %#ok<AGROW>
            wts(end+1) = t(end); %#ok<AGROW>
        end
    end
end
B = T(ismember(T.task, common), :);
M.logMeanNN  = log(mean(nn));
M.logRMSSD   = log(sqrt(mean(dd.^2)));
M.logMASD    = log(median(abs(dd)));   % outlier-resistant twin of RMSSD
M.pNN20      = 100 * mean(abs(dd) > 0.02);
M.rmssdRatio = sqrt(mean(dd.^2)) / max(median(abs(dd)), eps);  % QC: >3 means outlier-driven
M.logMedPV   = log(median(B.RR_max_min));
M.breathRate = 60 / median(B.len);
% individualized respiratory band
fR = 1 / median(B.len);
lo = max(0.04, 0.75*fR); hi = min(0.5, 1.25*fR);
hf = nan;
if ~isempty(wts)
    tot = 0; wsum = 0;
    for i = 1:numel(runs)
        v = runs(i).nn; f = v(isfinite(v));
        if numel(f) < 40, continue; end
        t = cumsum(f); tu = t(1):0.25:t(end);
        if numel(tu) < 128, continue; end
        x = detrend(interp1(t, f, tu, 'pchip'));
        w = min(256, 2^floor(log2(numel(x))));
        [px, fq] = pwelch(x, hamming(w), round(w/2), [], 4);
        sel = fq >= lo & fq < hi;
        if any(sel), tot = tot + trapz(fq(sel), px(sel)) * t(end); wsum = wsum + t(end); end
    end
    if wsum > 0, hf = tot / wsum; end
end
M.logHF = log(max(hf, eps));
% adjusted RSA level: fitted log PV at a fixed reference breath.
% bt(3) is the DEPTH slope - the transfer-function idea on the regressor that
% does not set the analysis window, so the order-statistics artifact does not
% apply (simulated null +0.07..0.10 on depth vs +0.49..0.82 on duration).
% bt(2), the duration coefficient, stays a nuisance term and is never a result.
y = log(B.RR_max_min); X = [ones(height(B),1), log(B.len), log(B.inhVol)];
if height(B) >= 20 && rank(X) == 3
    bt = X \ y;
    M.adjLogRSA  = bt(1) + bt(2)*log(LREF) + bt(3)*log(VREF);
    M.depthSlope = bt(3);
else
    M.adjLogRSA = NaN; M.depthSlope = NaN;
end

% evoked RSA: breath-phase-locked average RR waveform, peak-to-trough.
% Duration divides out by construction and unlocked low-frequency power
% averages toward zero as 1/sqrt(nBreaths), so there is no competing mechanism
% and the sign is unambiguous.
G = 100; W = nan(height(B), G);
for k = 1:height(B)
    a = B.finalOnset(k); z = min(a + round(B.len(k)*fs), numel(rrint));
    if z - a < 20, continue; end
    seg = rrint(a:z);
    W(k, :) = interp1(linspace(0, 1, numel(seg)), seg, linspace(0, 1, G));
end
mw = mean(W, 1, 'omitnan');
M.evokedRSA = max(mw) - min(mw);
end

function S = bootSE(runs, T, common, VREF, LREF, NBOOT, SEGSEC, METS, rrint, fs)
B = T(ismember(T.task, common), :);
vals = nan(NBOOT, numel(METS));
for b = 1:NBOOT
    rb = runs;
    for i = 1:numel(runs)
        v = runs(i).nn; f = find(isfinite(v));
        if numel(f) < 40, continue; end
        seg = max(10, round(SEGSEC / mean(v(f), 'omitnan')));
        nseg = max(1, floor(numel(v) / seg));
        pick = randi(nseg, 1, nseg);
        idx = [];
        for q = pick, idx = [idx, (q-1)*seg + (1:seg)]; end %#ok<AGROW>
        idx = idx(idx <= numel(v));
        rb(i).nn = v(idx);
    end
    Bb = B(randi(height(B), height(B), 1), :);
    try
        m = metrics(rb, Bb, common, VREF, LREF, rrint, fs);
        for k = 1:numel(METS), vals(b,k) = m.(METS{k}); end
    catch
    end
end
for k = 1:numel(METS), S.(METS{k}) = std(vals(:,k), 'omitnan'); end
end

function s = ternS(c, a, b), if c, s = a; else, s = b; end, end

function hb = readNPZ(f)
tmp = tempname; mkdir(tmp); c = onCleanup(@() rmdir(tmp, 's'));
unzip(f, tmp); d = dir(fullfile(tmp, 'hb.npy'));
fid = fopen(fullfile(d.folder, d.name), 'r');
fread(fid, 8, '*uint8'); hlen = fread(fid, 1, 'uint16'); fread(fid, hlen, '*char');
hb = fread(fid, inf, 'double'); fclose(fid);
end
