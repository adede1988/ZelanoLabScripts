function driver_rerun()
% DRIVER_RERUN  Re-run the breathing (audiobook+focusedBreathing) units of both
% the per-breath gamma and the 3-band spectrogram batches. The resume guards in
% run_gamma_all / run_spectro_all skip cue/thresh/O15 units whose outputs already
% exist, so only the deleted breathing outputs are regenerated (with the fixed
% focus-label mapping + p3 landmark). Reads finals from R: (mapped in the caller).
fprintf('driver_rerun start %s\n', datestr(now));
try, run_gamma_all(true);   catch e, fprintf('GAMMA FAIL: %s\n', e.message); end
try, run_spectro_all(true); catch e, fprintf('SPECTRO FAIL: %s\n', e.message); end
fprintf('DRIVER DONE %s\n', datestr(now));
end
