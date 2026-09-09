% chirp_smoke_battery2 -- clean full-battery discrimination check over 4 families (fast code).
setup_chirpAnalysis_paths(false);
C = chirp_config(); fs = C.fs;
C.ridge.nSurr = 40; C.ridge.nNullTrials = 2; C.phase.nPerm = 300;   % cheap-ish for the smoke
s0p = round((C.epochWin(1)-C.padSec)*fs); s1p = round((C.epochWin(2)+C.padSec)*fs);
off = s0p:s1p; coreIdx = off>=round(C.epochWin(1)*fs) & off<=round(C.epochWin(2)*fs); tMs = off(coreIdx)/fs*1000;

base = struct('overlapMs',800,'Df',7,'ampRatio',[1 1],'fLo',33,'snrDb',8,'burstMs',1500,'fastDropMs',400,'jitterMs',80,'seed',1);
fams = {'singleLin','singleNL','dualIndep','dualLocked'};
nTr = 14;
fprintf('\n%-11s | nR2 | phase permZ vtP | beat fSig mLSNR | chirp S/D/A | temp corr\n', 'family');
fprintf('%s\n', repmat('-',1,92));
for fi = 1:numel(fams)
    fam = fams{fi}; X = zeros(nTr, numel(off));
    for k = 1:nTr, P = base; P.seed = k*7 + fi*1000; [x,~] = chirp_synth_trial(fam, P, C); X(k,:) = x; end
    E = struct('dataPad',X,'coreIdx',coreIdx,'tMs',tMs,'valid',true(nTr,1),'fs',fs,'onsets',(1:nTr)','padSec',C.padSec,'epochWin',C.epochWin);
    tfr = chirp_tfr_faslt(E, C); R = chirp_ridge(tfr, E, C);
    ph = chirp_phase_test(E, R, C); be = chirp_beat_test(E, R, C); td = chirp_temporal_decoup(E, R, C);
    cls = cell(1,nTr);
    for i = 1:nTr
        wfit = R.trial(i).transitionWin; if isempty(wfit), wfit = R.trial(i).coexistWin; end
        o = chirp_chirplet(E.dataPad(i,coreIdx), tMs, wfit, C); cls{i} = o.classification;
    end
    nr2 = mean(arrayfun(@(t)t.nRidgesSig,R.trial)>1);
    fprintf('%-11s | %.2f | %6.2f %5.3f | %.2f %6.2f | %.2f/%.2f/%.2f | %5.2f (n=%d)\n', ...
        fam, nr2, ph.summary.perm_z, ph.summary.vtest_p_vs0, ...
        nanz(be.summary.fracSig), be.summary.medianLogSNR, ...
        mean(strcmp(cls,'single')), mean(strcmp(cls,'dual')), mean(strcmp(cls,'ambig')), ...
        td.summary.corrHiLo, td.summary.nUsed);
end
fprintf('\nExpect: single* -> high permZ, low beat fSig, chirp S high; dual* -> high nR2, beat fSig up, chirp D high.\n');
fprintf('BATTERY2_DONE\n');
function v = nanz(x), if isempty(x), v=NaN; else, v=x; end, end
