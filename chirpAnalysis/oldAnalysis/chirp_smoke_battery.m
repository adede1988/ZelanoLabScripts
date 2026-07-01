% chirp_smoke_battery -- run the FULL test battery on synthetic single vs dual; confirm it discriminates.
setup_chirpAnalysis_paths(false);
C = chirp_config(); fs = C.fs;
C.surr.n = 30; C.ridge.nSurr = 20; C.phase.nPerm = 200;   % cheap for the smoke

s0p = round((C.epochWin(1)-C.padSec)*fs); s1p = round((C.epochWin(2)+C.padSec)*fs);
off = s0p:s1p; coreIdx = off>=round(C.epochWin(1)*fs) & off<=round(C.epochWin(2)*fs);
tMs = off(coreIdx)/fs*1000;

base = struct('overlapMs',800,'Df',7,'ampRatio',[1 1],'fLo',33,'snrDb',8,'burstMs',1500,'fastDropMs',400,'jitterMs',60,'seed',1);
fams = {'singleLin','dualIndep'};
nTr = 10;
for fi = 1:numel(fams)
    fam = fams{fi};
    X = zeros(nTr, numel(off));
    for k = 1:nTr
        P = base; P.seed = k*7 + fi*1000;
        [x,~] = chirp_synth_trial(fam, P, C); X(k,:) = x;
    end
    E = struct('dataPad',X,'coreIdx',coreIdx,'tMs',tMs,'valid',true(nTr,1),'fs',fs, ...
        'onsets',(1:nTr)','padSec',C.padSec,'epochWin',C.epochWin);
    tfr = chirp_tfr_faslt(E, C);
    R   = chirp_ridge(tfr, E, C);
    ph  = chirp_phase_test(E, R, C);
    be  = chirp_beat_test(E, R, C);
    td  = chirp_temporal_decoup(E, R, C);
    cls = cell(1,nTr);
    for i = 1:nTr
        o = chirp_chirplet(E.dataPad(i,coreIdx), tMs, R.trial(i).transitionWin, C);
        cls{i} = o.classification;
    end
    nrs = arrayfun(@(t) t.nRidgesSig, R.trial);
    fprintf('\n=== %s (Df=%d overlap=%dms ampRatio=1:1 snr=%ddB) ===\n', fam, base.Df, base.overlapMs, base.snrDb);
    fprintf(' ridge: mean nRidgesSig=%.2f  frac>1=%.2f  anchors f_hi=%.1f f_lo=%.1f\n', ...
        mean(nrs,'omitnan'), mean(nrs>1), R.anchors.f_hi, R.anchors.f_lo);
    fprintf(' PHASE: perm_z=%.2f Rbar=%.2f vtest_p=%.3g nUsed=%d   (single->high perm_z)\n', ...
        ph.summary.perm_z, ph.summary.Rbar, ph.summary.vtest_p_vs0, ph.summary.nUsed);
    fprintf(' BEAT : medianLogSNR=%.2f fracSig=%.2f nUsed=%d fBeat=%.1f   (dual->high fracSig)\n', ...
        be.summary.medianLogSNR, be.summary.fracSig, be.summary.nUsed, be.params.fBeat);
    fprintf(' TEMP : corrHiLo=%.2f varLatDiff=%.0f nUsed=%d   (single->corr~1 low var)\n', ...
        td.summary.corrHiLo, td.summary.varLatDiff, td.summary.nUsed);
    fprintf(' CHIRP: fracSingle=%.2f fracDual=%.2f fracAmbig=%.2f\n', ...
        mean(strcmp(cls,'single')), mean(strcmp(cls,'dual')), mean(strcmp(cls,'ambig')));
end
fprintf('\nBATTERY_SMOKE_DONE\n');
