function rsa_robust()
% RMSSD squares successive differences, so a handful of large ones dominate it.
% pNN20 moving opposite to RMSSD is the signature of exactly that. Compare
% RMSSD against MASD (median absolute successive difference), which is
% outlier-resistant, on the four primary-contrast subjects.

P      = ohrv_config();

sp     = P.code;
outDir = P.work;
SPONT   = ["audio","focus","naturalFocus"]; DROPSEC = 30;
want = ["250929_Dupi_NMH_GH_1","250929_Dupi_NMH_GH_2","251002_Dupi_NMH_AB_1", ...
        "251002_Dupi_NMH_AB_2","251030_Dupi_NMH_DB_1","251030_Dupi_NMH_DB_2", ...
        "250818_Dupi_NMH_JH_1","250818_Dupi_NMH_JH_2"];

fprintf('%-24s %6s %8s %8s %8s %8s %8s\n', 'session','nDiff','RMSSD','MASD','p95|d|','max|d|','RMSSD/MASD');
R = table();
for w = want
    L = load(fullfile(outDir, w + "_slim.mat"));
    T = L.T; fs = L.fs; nSamp = L.nSamp; clear L
    hb = readNPZ(fullfile(outDir, w + "_beats.npz"));
    keep = T.goodBreath==1 & isfinite(T.RR_max_min) & T.RR_max_min>0 & T.len>=1.5 & T.len<=15 ...
         & isfinite(T.inhDur) & T.inhDur>0 & T.inhDur<T.len & isfinite(T.inhVol) & T.inhVol>0 ...
         & isfinite(T.finalOnset) & T.finalOnset>=1 & (T.finalOnset+round(T.len*fs))<=nSamp ...
         & strcmpi(T.noseMouth,"nose") & ismember(T.task, SPONT);
    T = T(keep,:);
    dd = [];
    for k = unique(T.condition)'
        m = T.condition==k;
        t0 = min(T.finalOnset(m)) + DROPSEC*fs; t1 = max(T.finalOnset(m)+round(T.len(m)*fs));
        bb = hb(hb>=t0 & hb<=t1); if numel(bb)<40, continue; end
        nn = diff(bb)/fs; ok = nn>0.3 & nn<2.0;
        md = movmedian(nn,11,'omitnan'); ok = ok & abs(nn-md)<=0.35*md;
        nn(~ok) = NaN; d = diff(nn); dd = [dd; d(isfinite(d))]; %#ok<AGROW>
    end
    rmssd = sqrt(mean(dd.^2)); masd = median(abs(dd));
    fprintf('%-24s %6d %8.4f %8.4f %8.4f %8.4f %8.2f\n', w, numel(dd), rmssd, masd, ...
        prctile(abs(dd),95), max(abs(dd)), rmssd/masd);
    R = [R; table(w, numel(dd), rmssd, masd, 'VariableNames', {'sess','n','RMSSD','MASD'})]; %#ok<AGROW>
end

fprintf('\n%-6s %12s %12s\n', 'subj','d_logRMSSD','d_logMASD');
for s = ["GH","AB","DB","JH"]
    i = find(contains(R.sess, "_"+s+"_"));
    fprintf('%-6s %+12.3f %+12.3f\n', s, ...
        log(R.RMSSD(i(2)))-log(R.RMSSD(i(1))), log(R.MASD(i(2)))-log(R.MASD(i(1))));
end
end

function hb = readNPZ(f)
tmp = tempname; mkdir(tmp); c = onCleanup(@() rmdir(tmp,'s'));
unzip(f,tmp); d = dir(fullfile(tmp,'hb.npy'));
fid = fopen(fullfile(d.folder,d.name),'r');
fread(fid,8,'*uint8'); hlen = fread(fid,1,'uint16'); fread(fid,hlen,'*char');
hb = fread(fid,inf,'double'); fclose(fid);
end
