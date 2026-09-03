function diag_excl()
% Why did 11 sessions collapse to a single block?
P      = ohrv_config();
outDir = P.work;
files = dir(fullfile(outDir, '*_slim.mat'));
fprintf('%-30s %6s %7s %7s %7s %7s %7s  %s\n', 'session','nRaw','good','nose','pass','nBlkRaw','nBlkOK','tasks(raw)');
for ii = 1:numel(files)
    L = load(fullfile(files(ii).folder, files(ii).name));
    T = L.T; fs = L.fs; nSamp = L.nSamp; sess = L.sessID; clear L
    nRaw = height(T);
    g    = sum(T.goodBreath == 1);
    nose = sum(strcmpi(T.noseMouth, "nose"));
    keep = T.goodBreath==1 & isfinite(T.RR_max_min) & T.RR_max_min>0 & T.len>=1.5 & T.len<=15 ...
         & isfinite(T.inhDur) & T.inhDur>0 & T.inhDur<T.len & isfinite(T.inhVol) ...
         & isfinite(T.finalOnset) & T.finalOnset>=1 & (T.finalOnset+round(T.len*fs))<=nSamp ...
         & strcmpi(T.noseMouth,"nose");
    K = T(keep,:);
    [ub,~,bi] = unique(K.condition); cnt = accumarray(bi,1);
    fprintf('%-30s %6d %7d %7d %7d %7d %7d  %s\n', sess, nRaw, g, nose, height(K), ...
        numel(unique(T.condition)), sum(cnt>=5), strjoin(unique(T.task)', ','));
end
end
