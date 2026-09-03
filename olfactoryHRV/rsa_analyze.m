function rsa_analyze()
% Per-session transfer-function gain, surrogate-calibrated, then change scores
% against an olfactory composite.

maxNumCompThreads(3);
rng(42);
P      = ohrv_config();
sp     = P.code;
outDir = P.work;
B      = 200;
MINBR  = 100;   % breaths per session
MINBLK = 5;     % breaths per block for that block to count

files = dir(fullfile(outDir, '*_slim.mat'));
fprintf('%d slim files\n\n', numel(files));

S = table();
tasksRows = table();

for ii = 1:numel(files)
    L = load(fullfile(files(ii).folder, files(ii).name));
    T = L.T; rrint = L.rrint; fs = L.fs; nSamp = L.nSamp; sessID = L.sessID; clear L
    bfile = fullfile(outDir, [sessID '_beats.npz']);
    if ~exist(bfile, 'file'), fprintf('%-32s no beats\n', sessID); continue; end
    hb = readNPZbeats(bfile);

    % ---------- exclusions ----------
    keep = T.goodBreath == 1 ...
         & isfinite(T.RR_max_min) & T.RR_max_min > 0 ...
         & T.len >= 1.5 & T.len <= 15 ...
         & isfinite(T.inhDur) & T.inhDur > 0 & T.inhDur < T.len ...
         & isfinite(T.inhVol) ...
         & isfinite(T.finalOnset) & T.finalOnset >= 1 ...
         & (T.finalOnset + round(T.len * fs)) <= nSamp ...
         & ~strcmpi(T.noseMouth, "mouth");
    % noseMouth is "nose", "NA" or "" (empty) in this cohort - no session carries
    % an explicit "mouth". Requiring =="nose" discarded 31-46% of breaths in 11
    % sessions, all in the audio/focus blocks, purely because the field was left
    % unpopulated. Respiration is recorded with a NASAL CANNULA, so a breath that
    % was detected at all with a valid RSA value was necessarily nasal; empty
    % means unrecorded metadata, not mouth breathing. The documented nose/mouth
    % manipulation occurs only in shadow blocks, which are paced and excluded
    % from the spontaneous panel regardless.
    T = T(keep, :);

    % blocks with too few breaths drop out (they cannot be centered usefully)
    [ub, ~, bi] = unique(T.condition);
    cnt = accumarray(bi, 1);
    good = ismember(T.condition, ub(cnt >= MINBLK));
    T = T(good, :); bi = bi(good);
    [~, ~, bi] = unique(bi);

    n = height(T); nBlk = numel(unique(bi));
    if n < MINBR || nBlk < 2
        fprintf('%-32s  EXCLUDED  n=%d nBlk=%d\n', sessID, n, nBlk);
        continue
    end

    % ---------- design ----------
    x1 = log(T.len);
    x2 = 1 - T.inhDur ./ T.len;
    v  = T.inhVol;
    x3 = (v - median(v)) ./ (1.4826 * mad(v, 1) + eps);
    X  = [blockCenter(x1, bi), blockCenter(x2, bi), blockCenter(x3, bi)];

    y  = blockCenter(log(T.RR_max_min), bi);
    bet = X \ y;

    % ---------- surrogate nulls ----------
    % NULL 1 (circular shift): keeps the RR series intact, only breaks its
    % alignment to the breaths. Vagal tone drifts slowly, so a shifted window
    % still contains real RSA at the subject's own breathing rate - this is a
    % CONSERVATIVE null that removes phase alignment, not respiratory
    % modulation.
    % NULL 2 (i.i.d. permuted NN): rebuilds an RR series from the same interval
    % distribution with all temporal structure destroyed. Its slope is the pure
    % order-statistics effect - longer window, more beats, wider range.
    % Their difference separates "mechanical" from "oscillation present but
    % misaligned".
    on = T.finalOnset; wl = round(T.len * fs);
    nullB = nan(B, 3); nullI = nan(B, 1);
    lo = round(30 * fs); hi = round(300 * fs);

    nn = diff(hb) / fs; nn = nn(nn > 0.3 & nn < 2.0);
    tg = (1:nSamp)' / fs;

    for s = 1:B
        rs  = circshift(rrint, randi([lo hi]));
        rmm = winRange(rs, on, wl, n);
        ok  = isfinite(rmm) & rmm > 0;
        if sum(ok) >= 0.9 * n
            nullB(s, :) = (X(ok, :) \ blockCenter(log(rmm(ok)), bi(ok)))';
        end

        if ~isempty(nn)
            p  = nn(randperm(numel(nn)));
            bt = cumsum(p);
            ri = interp1(bt, p, tg, 'linear', 'extrap');
            rmm2 = winRange(ri, on, wl, n);
            ok2  = isfinite(rmm2) & rmm2 > 0;
            if sum(ok2) >= 0.9 * n
                bb = X(ok2, :) \ blockCenter(log(rmm2(ok2)), bi(ok2));
                nullI(s) = bb(1);
            end
        end
    end
    nb     = nullB(isfinite(nullB(:, 1)), :);
    theta  = bet(1) - mean(nb(:, 1));      % DURATION slope, surrogate-corrected
    seTh   = std(nb(:, 1));
    nullIm = mean(nullI, 'omitnan');
    % DEPTH slope. Breath depth does not set the length of the analysis window,
    % so the order-statistics artifact that contaminates the duration slope does
    % not apply here - this is the better-posed version of the same question,
    % "how much cardiac modulation does deeper breathing recruit".
    thetaV = bet(3) - mean(nb(:, 3));
    seThV  = std(nb(:, 3));
    thetaE = bet(2) - mean(nb(:, 2));
    seThE  = std(nb(:, 2));

    midRR = (T.maxRR + T.minRR) / 2;
    tok   = split(string(sessID), '_');
    subj  = tok(4); if subj == "TPB", subj = "TB"; end

    S = [S; table(string(sessID), subj, double(tok(5)), n, nBlk, ...
        bet(1), bet(2), bet(3), mean(nb(:,1)), nullIm, mean(nb(:,3)), ...
        theta, seTh, thetaV, seThV, thetaE, seThE, ...
        60/median(midRR), std(x1), median(T.RR_max_min), ...
        'VariableNames', {'sessID','subj','sess','n','nBlk', ...
        'b_len','b_exh','b_vol','nullShift','nullIID','nullVol', ...
        'theta','seTheta','thetaVol','seThetaVol','thetaExh','seThetaExh', ...
        'HR','sdLogLen','medRSA'})]; %#ok<AGROW>

    fprintf(['%-30s n=%4d blk=%2d | dur b=%+.3f nullS=%+.3f nullI=%+.3f th=%+.3f(%.3f)' ...
             ' | depth b=%+.3f null=%+.3f th=%+.3f(%.3f) | HR=%4.1f\n'], ...
        sessID, n, nBlk, bet(1), mean(nb(:,1)), nullIm, theta, seTh, ...
        bet(3), mean(nb(:,3)), thetaV, seThV, 60/median(midRR));

    % ---------- length x condition diagnostic ----------
    ut = unique(T.task);
    for k = 1:numel(ut)
        m = T.task == ut(k);
        if sum(m) >= 30
            [~, ~, bsub] = unique(bi(m));           % re-index so accumarray has no empty bins
            bb = blockCenter(x1(m), bsub);
            yy = blockCenter(log(T.RR_max_min(m)), bsub);
            if std(bb) > 0
                tasksRows = [tasksRows; table(string(sessID), subj, ut(k), sum(m), (bb\yy), ...
                    'VariableNames', {'sessID','subj','task','n','b_len'})]; %#ok<AGROW>
            end
        end
    end
    clear T rrint
end

writetable(S, fullfile(outDir, 'session_gains.csv'));
writetable(tasksRows, fullfile(outDir, 'task_slopes.csv'));

fprintf('\nANALYZE DONE  (%d sessions)\n', height(S));
end


function xc = blockCenter(x, bi)
mu = accumarray(bi, x, [], @mean);
xc = x - mu(bi);
end

function rmm = winRange(r, on, wl, n)
rmm = nan(n, 1);
L = numel(r);
for k = 1:n
    a = on(k); b = min(on(k) + wl(k), L);
    if a >= 1 && b > a
        w = r(a:b);
        rmm(k) = max(w) - min(w);
    end
end
end

function hb = readNPZbeats(f)
tmp = tempname; mkdir(tmp); c = onCleanup(@() rmdir(tmp, 's'));
unzip(f, tmp); d = dir(fullfile(tmp, 'hb.npy'));
fid = fopen(fullfile(d.folder, d.name), 'r');
fread(fid, 8, '*uint8'); hlen = fread(fid, 1, 'uint16'); fread(fid, hlen, '*char');
hb = fread(fid, inf, 'double'); fclose(fid);
end
