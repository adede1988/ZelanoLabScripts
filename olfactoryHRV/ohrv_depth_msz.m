function M = ohrv_depth_msz()
%OHRV_DEPTH_MSZ  Per-session depth slope in interpretable units: ms per z.
%
%   Same design as rsa_analyze's per-session fit -
%       outcome ~ log(breath length) + exhale dominance + robustZ(inhaled volume)
%   block-centred within session - but with a LINEAR outcome in milliseconds
%   instead of log(RR_max_min). The depth coefficient is then literally
%   "milliseconds of respHRV per standard deviation of breath depth", which is
%   what a reader expects from an axis labelled ms/z.
%
%   rsa_analyze's b_vol uses a log outcome, so its units are log-RSA per z, NOT
%   ms/z. The two change scores agree closely (r = 0.95 across the primary
%   contrast), so this is a relabelling of the same effect in readable units
%   rather than a different result.
%
%   Returns a table: sessID, subj, sess, msz.

P = ohrv_config(); outDir = P.work;
G = readtable(fullfile(outDir,'report_sessions.csv'),'TextType','string');

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
    y  = bc(T.RR_max_min*1000, bi);                       % linear outcome, ms
    bet = X \ y;
    msz(ii) = bet(3);                                     % ms per z of depth
end

M = table(G.sessID, G.subj, G.sess, msz, ...
          'VariableNames', {'sessID','subj','sess','msz'});
end

function xc = bc(x, bi)
mu = accumarray(bi, x, [], @mean); xc = x - mu(bi);
end
