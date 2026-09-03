function rsa_secheck()
% Cross-check the surrogate-derived SE against OLS and block-clustered SEs,
% then recompute kappa under each. kappa is the whole conclusion, so the SE
% estimate it rests on needs to be verified rather than assumed.

maxNumCompThreads(3);
P      = ohrv_config();
sp     = P.code;
outDir = P.work;
MINBR = 100; MINBLK = 5;

G = readtable(fullfile(outDir, 'session_gains.csv'), 'TextType', 'string');
R = table();

for ii = 1:height(G)
    L = load(fullfile(outDir, G.sessID(ii) + "_slim.mat"));
    T = L.T; fs = L.fs; nSamp = L.nSamp; clear L

    keep = T.goodBreath == 1 & isfinite(T.RR_max_min) & T.RR_max_min > 0 ...
         & T.len >= 1.5 & T.len <= 15 ...
         & isfinite(T.inhDur) & T.inhDur > 0 & T.inhDur < T.len & isfinite(T.inhVol) ...
         & isfinite(T.finalOnset) & T.finalOnset >= 1 ...
         & (T.finalOnset + round(T.len*fs)) <= nSamp & strcmpi(T.noseMouth, "nose");
    T = T(keep, :);
    [ub, ~, bi] = unique(T.condition); cnt = accumarray(bi, 1);
    T = T(ismember(T.condition, ub(cnt >= MINBLK)), :);
    [~, ~, bi] = unique(T.condition);
    n = height(T); if n < MINBR || numel(unique(bi)) < 2, continue; end

    x1 = log(T.len); x2 = 1 - T.inhDur./T.len;
    v = T.inhVol; x3 = (v - median(v)) ./ (1.4826*mad(v,1) + eps);
    X = [bc(x1,bi), bc(x2,bi), bc(x3,bi)];
    y = bc(log(T.RR_max_min), bi);

    bet = X \ y; e = y - X*bet;
    XtXi = inv(X' * X);
    df   = n - numel(unique(bi)) - 3;
    seOLS = sqrt(XtXi(1,1) * (e'*e) / df);

    % block-clustered sandwich
    M = zeros(3);
    for g = unique(bi)'
        m = bi == g; Xg = X(m,:); eg = e(m);
        M = M + (Xg'*eg) * (Xg'*eg)';
    end
    nc = numel(unique(bi));
    adj = (nc/(nc-1)) * ((n-1)/df);
    Vc  = XtXi * M * XtXi * adj;
    seCL = sqrt(Vc(1,1));

    R = [R; table(G.sessID(ii), G.subj(ii), G.sess(ii), n, numel(unique(bi)), ...
        G.theta(ii), G.seTheta(ii), seOLS, seCL, ...
        'VariableNames', {'sessID','subj','sess','n','nBlk','theta','seNull','seOLS','seCluster'})]; %#ok<AGROW>
    clear T
end

disp(R)
fprintf('\nmedian seNull=%.3f  seOLS=%.3f  seCluster=%.3f\n', ...
    median(R.seNull), median(R.seOLS), median(R.seCluster));

% recompute kappa under each SE choice
D = readtable(fullfile(outDir, 'change_scores.csv'), 'TextType', 'string');
names = {'seNull','seOLS','seCluster'};
fprintf('\nkappa = mean SE^2(dTheta) / var(dTheta)   [var = %.4f over %d intervals]\n', ...
    var(D.dTheta), height(D));
for k = 1:3
    sd = nan(height(D),1);
    for j = 1:height(D)
        a = R.subj == D.subj(j) & R.sess == D.from(j);
        b = R.subj == D.subj(j) & R.sess == D.to(j);
        sd(j) = sqrt(R.(names{k})(a)^2 + R.(names{k})(b)^2);
    end
    kap = mean(sd.^2) / var(D.dTheta);
    fprintf('  %-10s  mean SE(dTheta)=%.3f   kappa=%.3f   attenuation=%.3f\n', ...
        names{k}, mean(sd), kap, 1/sqrt(1+kap));
end
writetable(R, fullfile(outDir, 'se_comparison.csv'));
end

function xc = bc(x, bi)
mu = accumarray(bi, x, [], @mean); xc = x - mu(bi);
end
