function sub = chirp_phase_test(E, R, C)
% CHIRP_PHASE_TEST  Segmented phase-continuity test (spec 6.5) -- PRIMARY single-channel test.
%   sub = chirp_phase_test(E, R, C)
%   E : epoch struct (chirp_epoch) for the OB channel (bestMac); E.dataPad/coreIdx/tMs/valid
%   R : ridge result (chirp_ridge); per-trial R.trial(i).fhat (POWER ridge, on R.tMs) + .transitionWin
%   C : chirp_config
%
%   For one coherent oscillator the broadband Hilbert phase increment matches the ridge-IF
%   integral in EVERY short sub-window -> residual concentrated near 0 (strong-for-SINGLE).
%   Two generators inject independent phase -> dispersed residual (inconclusive). f-hat is the
%   POWER ridge, NEVER dphi/dt (rule a, asserted). Long-burst correction: segment the transition
%   into ~150 ms sub-windows so integrated IF-bias << 1 cycle.
%
%   sub.params / sub.trial(i) / sub.summary  (spec 7.1)

    assert(strcmp(C.phase.ifSource,'ridgePower'), ...
        'chirp_phase_test:ifSource', 'f-hat MUST be the power ridge (rule a), got %s', C.phase.ifSource);
    fs = E.fs; tMs = E.tMs;
    assert(isequal(tMs, R.tMs), 'chirp_phase_test:grid', 'epoch/ridge time grids differ');

    % broadband-gamma analytic phase on the PADDED epoch, then trim (unwrap before trim)
    bb = chirp_bbfilt(E.dataPad, fs, C.bbBand, C);
    bb = bb - mean(bb, 2, 'omitnan');
    z = hilbert(bb')';                 % analytic per row
    phiPad = unwrap(angle(z), [], 2);
    envPad = abs(z);
    phi = phiPad(:, E.coreIdx);
    env = envPad(:, E.coreIdx);

    subWin = C.phase.subWinMs;
    nTrial = size(E.dataPad,1);
    trial = struct('subWinResid',{},'meanRbar',{},'slipLatencies',{},'included',{},'reason',{});
    residByTrial = {};

    % baseline power for the existence gate (pre-sniff [-500 -100] ms)
    blMask = tMs >= -500 & tMs <= -100;

    for i = 1:nTrial
        T = struct('subWinResid',[],'meanRbar',NaN,'slipLatencies',[],'included',false,'reason','');
        if ~E.valid(i) || isempty(R.trial(i).fhat), T.reason='noTrialOrRidge'; trial(i)=T; continue; end
        tw = R.trial(i).transitionWin;
        if isempty(tw), tw = R.trial(i).coexistWin; end    % true-dual has no sweep -> use coexistence
        if isempty(tw) || (tw(2)-tw(1)) < C.phase.minTransSec*1000, T.reason='shortWin'; trial(i)=T; continue; end

        fhat = R.trial(i).fhat;
        % existence gate (rule b: response existence only, never shape)
        if C.phase.requireResponse
            inTr = tMs>=tw(1) & tMs<=tw(2);
            blP = median(env(i,blMask).^2,'omitnan'); trP = median(env(i,inTr).^2,'omitnan');
            if ~(trP > C.phase.minRespZ * blP), T.reason='noResponse'; trial(i)=T; continue; end
        end
        % drop if ridge NaN over too much of the transition
        inTr = tMs>=tw(1) & tMs<=tw(2);
        if mean(isnan(fhat(inTr))) > C.phase.maxRidgeGap, T.reason='ridgeGap'; trial(i)=T; continue; end

        % tile transition into consecutive sub-windows
        edges = tw(1):subWin:tw(2);
        if numel(edges) < 3, T.reason='shortTransition'; trial(i)=T; continue; end
        nsub = numel(edges)-1;
        resid = nan(1,nsub); subCtr = nan(1,nsub);
        for k = 1:nsub
            a = edges(k); b = edges(k+1); subCtr(k) = (a+b)/2;
            ia = find(tMs>=a,1,'first'); ib = find(tMs<=b,1,'last');
            if isempty(ia)||isempty(ib)||ib<=ia, continue; end
            phiObsInc = phi(i,ib) - phi(i,ia);
            fseg = fhat(ia:ib); tseg = tMs(ia:ib)/1000;
            if any(isnan(fseg)), continue; end
            phiPredInc = 2*pi*trapz(tseg, fseg);
            resid(k) = wrapToPi(phiObsInc - phiPredInc);
        end
        good = ~isnan(resid);
        if nnz(good) < 2, T.reason='fewSubWin'; trial(i)=T; continue; end
        T.subWinResid = resid(good);
        T.meanRbar = abs(mean(exp(1i*resid(good))));
        T.slipLatencies = subCtr(good & abs(resid) > C.phase.slipThreshRad);
        T.included = true;
        trial(i) = T;
        residByTrial{end+1} = resid(good); %#ok<AGROW>
    end

    % pool across trials
    Rall = [residByTrial{:}];
    summary = struct('nUsed',numel(residByTrial),'Rbar',NaN,'rayleigh_z',NaN,'rayleigh_p',NaN, ...
        'vtest_p_vs0',NaN,'perm_z',NaN,'perm_p',NaN,'perm_z0',NaN,'perm_p0',NaN,'meanResidDeg',NaN, ...
        'crossTrialR',NaN,'crossTrial_p',NaN,'crossTrial_n',0);
    if ~isempty(Rall)
        n = numel(Rall); Rbar = abs(mean(exp(1i*Rall)));
        summary.Rbar = Rbar; summary.rayleigh_z = n*Rbar^2;
        summary.rayleigh_p = exp(-summary.rayleigh_z)*(1 + (2*summary.rayleigh_z-summary.rayleigh_z^2)/(4*n));
        summary.meanResidDeg = rad2deg(angle(mean(exp(1i*Rall))));
        % V-test against preferred direction 0
        V = sum(cos(Rall)); Rbar0 = V/n; u = Rbar0*sqrt(2*n);
        summary.vtest_p_vs0 = 1 - normcdf(u);
        % per-session rotation permutation null (absorbs n-dependence).
        %  perm_z  = concentration ANYWHERE (Rbar);  perm_z0 = concentration AT ZERO (mean cos resid,
        %  the SINGLE-specific signature). A true dual concentrates at a NONZERO offset -> low perm_z0.
        R0obs = mean(cos(Rall));
        nP = C.phase.nPerm; rng(C.rngSeed,'twister');
        nullR = zeros(1,nP); nullR0 = zeros(1,nP);
        for b = 1:nP
            rr = [];
            for c = 1:numel(residByTrial)
                rr = [rr, wrapToPi(residByTrial{c} + 2*pi*rand)]; %#ok<AGROW>
            end
            nullR(b) = abs(mean(exp(1i*rr)));
            nullR0(b) = mean(cos(rr));
        end
        summary.perm_z = (Rbar - mean(nullR))/max(std(nullR),eps);
        summary.perm_p = (1 + sum(nullR >= Rbar))/(1 + nP);
        summary.perm_z0 = (R0obs - mean(nullR0))/max(std(nullR0),eps);
        summary.perm_p0 = (1 + sum(nullR0 >= R0obs))/(1 + nP);

        % CROSS-TRIAL consistency (the spec's actual logic, robust to a constant ridge-IF bias):
        % one oscillator -> each trial's mean residual is the SAME (clustered across trials) -> high
        % crossTrialR; independent-phase duals scatter per-trial means -> low. This is the
        % discriminating single-confirmer (concentration AT ZERO fails under ridge bias).
        ptm = cellfun(@(r) angle(mean(exp(1i*r))), residByTrial);
        nTrR = numel(ptm);
        if nTrR >= 3
            crossR = abs(mean(exp(1i*ptm)));
            summary.crossTrialR = crossR; summary.crossTrial_n = nTrR;
            zc = nTrR*crossR^2;
            summary.crossTrial_p = exp(-zc)*(1 + (2*zc - zc^2)/(4*nTrR));
        else
            summary.crossTrialR = NaN; summary.crossTrial_n = nTrR; summary.crossTrial_p = NaN;
        end
    end

    sub.params = struct('f_hi',R.anchors.f_hi,'f_lo',R.anchors.f_lo,'subWinLen',subWin, ...
        'filterSpec',C.bbBand,'ifSource',C.phase.ifSource);
    sub.trial = trial;
    sub.summary = summary;
end
