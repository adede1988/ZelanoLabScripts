% chirp_dbg_chirplet2 -- fast classifier check over 4 families (oracle window, no ridge/surrogates).
setup_chirpAnalysis_paths(false);
C = chirp_config(); fs = C.fs;
s0p = round((C.epochWin(1)-C.padSec)*fs); s1p = round((C.epochWin(2)+C.padSec)*fs);
off = s0p:s1p; coreIdx = off>=round(C.epochWin(1)*fs) & off<=round(C.epochWin(2)*fs); tMs = off(coreIdx)/fs*1000;
fams = {'singleLin','singleNL','dualIndep','dualLocked'};
N = 12;
fprintf('\n%-11s  fracSingle fracDual fracAmbig  (oracle coexistWin, snr=10)\n', 'family');
for fi = 1:numel(fams)
    cls = cell(1,N);
    for k = 1:N
        P = struct('overlapMs',800,'Df',7,'ampRatio',[1 1],'fLo',33,'snrDb',10,'burstMs',1500,'fastDropMs',400,'jitterMs',60,'seed',k*7+fi*100);
        [x,gt] = chirp_synth_trial(fams{fi}, P, C);
        w = gt.coexistWinMs;                       % pipeline mostly uses coexistWin (ridge rarely finds a sweep)
        o = chirp_chirplet(x(coreIdx), tMs, w, C);
        cls{k} = o.classification;
    end
    fprintf('%-11s   %.2f      %.2f     %.2f\n', fams{fi}, ...
        mean(strcmp(cls,'single')), mean(strcmp(cls,'dual')), mean(strcmp(cls,'ambig')));
end
fprintf('DBG2_DONE\n');
