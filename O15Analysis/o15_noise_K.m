function K = o15_noise_K()
% O15_NOISE_K  Single source of truth for the RELATIVE sharp-deflection threshold
%   K (in robust-SDs) used by the O15 noise rejection. This is the cue analysis's
%   relative rule (cue_noise_trials, robust-z zd>K) ported to O15: a 10 ms window
%   is an "event" if its max-min range sits >K robust-SDs above that channel's
%   TYPICAL 10 ms swing, and a sniff is rejected if its epoch contains any event.
%
%   CALIBRATED for the O15 dataset (o15_calibrate_noise, 2026-06-29, 2041 sniffs /
%   26 sessions): K=15 drops ~16.5% of sniffs dataset-wide with NO session fully
%   excluded (worst ~95%). (cue's K=10 drops only 8% on cue but ~31.6% on O15 -
%   O15 bestMac channels swing much harder - and fully excludes KS_3, so it is NOT
%   reused here.) Other points: K=20 -> 9.8%, K=12 -> 22.9%, K=10 -> 31.6%.
%
%   Both run_o15_ztfr and run_o15_gamma_epochs read K from here, so changing the
%   threshold is a one-line edit and stays consistent across the pipeline.

    K = 15;
end
