function check_msz()
%CHECK_MSZ  Is the depth slope in ms/z equivalent to the log-outcome version?
%
%   FIG2 currently plots d_b_vol, the depth coefficient from
%       log(RR_max_min) ~ log(len) + exhaleDom + robustZ(inhVol)
%   whose units are log-RSA per z, NOT ms/z. This recomputes the same design
%   with a LINEAR outcome in ms so the coefficient really is ms per z, and
%   compares the resulting change-vs-olfaction correlation.

P = ohrv_config(); outDir = P.work;
G = readtable(fullfile(outDir,'report_sessions.csv'),'TextType','string');
D = readtable(fullfile(outDir,'report_changes.csv'),'TextType','string');
Q = D(D.from==1 & D.to==2,:);

msz = nan(height(G),1);
for ii = 1:height(G)
    f = fullfile(outDir, G.sessID(ii) + "_slim.mat");
    if ~exist(f,'file'), continue; end
    L = load(f); T = L.T; fs = L.fs; nSamp = L.nSamp; clear L
    keep = T.goodBreath==1 & isfinite(T.RR_max_min) & T.RR_max_min>0 ...
         & T.len>=1.5 & T.len<=15 & isfinite(T.inhDur) & T.inhDur>0 & T.inhDur<T.len ...
         & isfinite(T.inhVol) & T.inhVol>0 & isfinite(T.finalOnset) & T.finalOnset>=1 ...
         & (T.finalOnset+round(T.len*fs))<=nSamp & ~strcmpi(T.noseMouth,"mouth");
    T = T(keep,:);
    [ub,~,bi] = unique(T.condition); cnt = accumarray(bi,1);
    T = T(ismember(T.condition, ub(cnt>=5)), :);
    [~,~,bi] = unique(T.condition);
    if height(T) < 100 || numel(unique(bi)) < 2, continue; end

    x1 = log(T.len);
    x2 = 1 - T.inhDur ./ T.len;
    v  = T.inhVol;
    x3 = (v - median(v)) ./ (1.4826*mad(v,1) + eps);      % robust z of depth
    X  = [bc(x1,bi) bc(x2,bi) bc(x3,bi)];
    y  = bc(T.RR_max_min*1000, bi);                       % LINEAR outcome, ms
    bet = X \ y;
    msz(ii) = bet(3);                                     % ms per z of depth
end
G.msz = msz;

fprintf('%-26s %6s %10s %10s\n','session','sess','b_vol(log)','msz(ms/z)');
for ii = 1:height(G)
    fprintf('%-26s %6d %10.3f %10.2f\n', G.sessID(ii), G.sess(ii), G.b_vol(ii), G.msz(ii));
end

% change scores over the primary contrast
dm = nan(height(Q),1);
for i = 1:height(Q)
    a = G.msz(G.subj==Q.subj(i) & G.sess==Q.from(i));
    b = G.msz(G.subj==Q.subj(i) & G.sess==Q.to(i));
    if ~isempty(a) && ~isempty(b), dm(i) = b - a; end
end
Q.d_msz = dm;

g = isfinite(Q.d_msz) & isfinite(Q.dOlf);
fprintf('\n--- session 1 -> 2, n = %d ---\n', sum(g));
fprintf('%-6s %12s %12s %10s\n','subj','d_b_vol(log)','d_msz(ms/z)','dOlf');
for i = find(g)'
    fprintf('%-6s %12.4f %12.2f %10.3f\n', Q.subj(i), Q.d_b_vol(i), Q.d_msz(i), Q.dOlf(i));
end
fprintf('\n  log version : Pearson %+.3f   Spearman %+.3f\n', ...
    corr(Q.dOlf(g),Q.d_b_vol(g)), corr(Q.dOlf(g),Q.d_b_vol(g),'Type','Spearman'));
fprintf('  ms/z version: Pearson %+.3f   Spearman %+.3f\n', ...
    corr(Q.dOlf(g),Q.d_msz(g)),   corr(Q.dOlf(g),Q.d_msz(g),'Type','Spearman'));
fprintf('  agreement between the two Deltas: r = %+.3f\n', corr(Q.d_b_vol(g),Q.d_msz(g)));

writetable(Q, fullfile(outDir,'msz_compare.csv'));
end

function xc = bc(x, bi)
mu = accumarray(bi, x, [], @mean); xc = x - mu(bi);
end
