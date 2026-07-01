function sub = chirp_spatial(macSig, fs, onsets, goodMac, R, C)
% CHIRP_SPATIAL  Hi-vs-lo spatial-profile similarity across macBP contacts (spec 6.8).
%   sub = chirp_spatial(macSig, fs, onsets, goodMac, R, C)
%   macSig  : [nMac x nSamp] continuous, ALL macBP* rows (NOT scalp/EEG)
%   goodMac : [nMac x 1] logical QC mask (channels passing upstream QC)
%   R       : ridge result (per-trial transitionWin from bestMac; anchors f_hi/f_lo)
%   C       : chirp_config
%
%   Per user D2: run across ALL good macBP, NO geometry/spacing gating (eligible iff >=2 good).
%   Per trial & contact: early-high-band power (f_hi, early window) and late-low-band power
%   (f_lo, late window); profileSim = corr of the hi vs lo spatial vectors across contacts. Two
%   generators -> low similarity; one oscillator -> identical early/late topography (high sim).

    fHi = R.anchors.f_hi; fLo = R.anchors.f_lo; bh = C.spatial.bandHalfHz;
    gi = find(goodMac(:)'); nG = numel(gi);
    nTrial = numel(onsets);

    sub.params = struct('fHiBand',[fHi-bh fHi+bh],'fLoBand',[fLo-bh fLo+bh], ...
        'earlyFrac',[0 1/3],'lateFrac',[2/3 1]);
    sub.trial = struct('profileHi',{},'profileLo',{},'profileSim',{},'included',{});
    eligible = nG >= C.spatial.minGoodContacts && isfinite(fHi) && isfinite(fLo) && abs(fHi-fLo)>=1;
    sub.summary = struct('eligible',eligible,'nGoodContacts',nG,'contactSpacing',NaN, ...
        'meanProfileSim',NaN,'nUsed',0);
    if ~eligible
        sub.trial(nTrial) = struct('profileHi',[],'profileLo',[],'profileSim',NaN,'included',false);
        return;
    end

    % build early-high and late-low band power per good contact, epoch-locked
    hiPow = cell(1,nG); loPow = cell(1,nG); Es = cell(1,nG);
    for jj = 1:nG
        m = gi(jj);
        Em = chirp_epoch(macSig(m,:), fs, onsets, C); Es{jj} = Em;
        hP = abs(hilbert(chirp_bbfilt(Em.dataPad, fs, clampB([fHi-bh fHi+bh],fs), C)')').^2;
        lP = abs(hilbert(chirp_bbfilt(Em.dataPad, fs, clampB([fLo-bh fLo+bh],fs), C)')').^2;
        hiPow{jj} = hP(:, Em.coreIdx); loPow{jj} = lP(:, Em.coreIdx);
    end
    tMs = Es{1}.tMs;

    sims = [];
    for i = 1:nTrial
        T = struct('profileHi',[],'profileLo',[],'profileSim',NaN,'included',false);
        cw = R.trial(i).coexistWin; tw = R.trial(i).transitionWin;
        w = tw; if isempty(w), w = cw; end
        if isempty(w) || ~Es{1}.valid(i), sub.trial(i)=T; continue; end
        eW = [w(1), w(1)+(w(2)-w(1))/3];           % early third
        lW = [w(1)+2*(w(2)-w(1))/3, w(2)];         % late third
        eM = tMs>=eW(1)&tMs<=eW(2); lM = tMs>=lW(1)&tMs<=lW(2);
        pHi = nan(nG,1); pLo = nan(nG,1);
        for jj = 1:nG
            if ~Es{jj}.valid(i), continue; end
            pHi(jj) = mean(hiPow{jj}(i,eM),'omitnan');
            pLo(jj) = mean(loPow{jj}(i,lM),'omitnan');
        end
        ok = isfinite(pHi)&isfinite(pLo);
        if nnz(ok) < C.spatial.minGoodContacts, sub.trial(i)=T; continue; end
        T.profileHi = pHi(ok)'; T.profileLo = pLo(ok)';
        if nnz(ok) >= 2
            cc = corrcoef(pHi(ok), pLo(ok)); T.profileSim = cc(1,2);
        end
        T.included = isfinite(T.profileSim);
        sub.trial(i) = T;
        if T.included, sims(end+1) = T.profileSim; end %#ok<AGROW>
    end
    sub.summary.meanProfileSim = meanSafe(sims);
    sub.summary.nUsed = numel(sims);
end

function b = clampB(b, fs), b(1)=max(b(1),1); b(2)=min(b(2),fs/2-1); end
function v = meanSafe(x), if isempty(x), v=NaN; else, v=mean(x); end, end
