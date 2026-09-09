% chirp_dbg_chirplet -- inspect chirplet atoms for one single + one dual trial (fixed mode, fast).
setup_chirpAnalysis_paths(false);
C = chirp_config(); fs = C.fs; C.ridge.nSurr = 10;   % cheap ridge for debug
s0p = round((C.epochWin(1)-C.padSec)*fs); s1p = round((C.epochWin(2)+C.padSec)*fs);
off = s0p:s1p; coreIdx = off>=round(C.epochWin(1)*fs) & off<=round(C.epochWin(2)*fs); tMs = off(coreIdx)/fs*1000;

for fam = {'singleLin','dualIndep'}
    P = struct('overlapMs',800,'Df',7,'ampRatio',[1 1],'fLo',33,'snrDb',12,'burstMs',1500,'fastDropMs',400,'jitterMs',0,'seed',5);
    [x,gt] = chirp_synth_trial(fam{1}, P, C);
    E = struct('dataPad',x,'coreIdx',coreIdx,'tMs',tMs,'valid',true,'fs',fs,'onsets',1,'padSec',C.padSec,'epochWin',C.epochWin);
    tfr = chirp_tfr_faslt(E, C); R = chirp_ridge(tfr, E, C);
    tw = R.trial(1).transitionWin; cw = R.trial(1).coexistWin;
    fprintf('\n=== %s ===  gt.transWin=[%.0f %.0f]  ridge transWin=[%s]  coexist=[%s]\n', fam{1}, ...
        gt.transitionWinMs(1), gt.transitionWinMs(2), num2str(round(tw)), num2str(round(cw)));
    % run chirplet with surrogate mode to also see the threshold, but cap nSurr
    Cd = C; Cd.chirplet.sigMode='surrogate'; Cd.chirplet.nSurr=30;
    wfit = tw; if isempty(wfit), wfit = cw; end
    o = chirp_chirplet(x(coreIdx), tMs, wfit, Cd);
    fprintf(' nAtomsSig=%d traj=%.2f class=%s c1=%.2f atomThr=%.3f reason=%s\n', ...
        o.nAtomsSig, o.trajConnectivity, o.classification, o.c1, o.atomThr, o.reason);
    if ~isempty(o.atoms)
        fprintf('   tc_s    fc_Hz  chirpHz/s  sigma_s   amp     eFrac\n');
        for k=1:size(o.atoms,1)
            fprintf('   %6.3f  %5.1f  %8.2f  %7.3f  %6.2f  %5.3f\n', o.atoms(k,1),o.atoms(k,2),o.atoms(k,3),o.atoms(k,4),o.atoms(k,5),o.atoms(k,7));
        end
    end
    % also fixed mode
    of = chirp_chirplet(x(coreIdx), tMs, wfit, C);
    fprintf(' [fixed thr=%.2f] nAtomsSig=%d class=%s\n', C.chirplet.atomEnergyThr, of.nAtomsSig, of.classification);
end
fprintf('DBG_DONE\n');
