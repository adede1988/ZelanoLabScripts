function run_cue_singletrial(sessFilter)
% RUN_CUE_SINGLETRIAL  (analysis2/3) Standalone re-plot of the single-trial raw
%   bestMac noise figure with sharp-deflection trials flagged red. The full
%   analysis3 run (run_cue_ztfr) also produces this figure; this is a thin
%   wrapper for re-plotting alone. Uses the shared cue_noise_trials /
%   cue_plot_singletrial so the noise labels match the spectrogram rejection.

    if nargin < 1, sessFilter = []; end
    cue_init_paths(); L = labPaths(); sep = 50;

    T = cue_session_table(false); T = T(T.onDisk, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n== %d/%d %s ==\n', i, height(T), id);
        try
            Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
            if ~isfield(od,'bestMac') || isempty(od.bestMac), fprintf('  no bestMac -> skip\n'); clear od; continue; end
            labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
            ci = find(strcmp(od.bestMac, labs), 1); if isempty(ci), clear od; continue; end
            sig = double(od.data(ci, :)); fs = od.fs;

            NT = cue_noise_trials(sig, fs, od.TTL.trialStart);
            figDir = fullfile(L.figPath, id, 'cueTask'); if ~isfolder(figDir), mkdir(figDir); end
            cue_plot_singletrial(NT, id, od.bestMac, sep, fullfile(figDir, 'singleTrialRawMac.png'));
            fprintf('  saved singleTrialRawMac.png (n=%d, %d noise)\n', sum(NT.ok), sum(NT.noisy(NT.ok)));
            clear od sig NT;
        catch ME
            fprintf('  FAILED: %s\n', ME.message); try, clear od; catch, end
        end
    end
    fprintf('\nsingle-trial raw pass done.\n');
end
