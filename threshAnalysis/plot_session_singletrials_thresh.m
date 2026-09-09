function plot_session_singletrials_thresh(sessID, outSub, epWin)
% PLOT_SESSION_SINGLETRIALS_THRESH  One figure PER TRIAL of a session's bipolar
%   (macBP) channels, TTL.start-locked, stacked. bestMac is drawn slightly bolder
%   and the RELATIVE sharp-deflection windows (zd>K, the same rule
%   thresh_noise_trials uses for rejection) are overlaid in red on bestMac. The
%   title says whether the whole trial was rejected.
%
%   plot_session_singletrials_thresh(sessID [, outSub, epWin])
%   Saves <figs>\<id>\threshTask\<outSub>\trial_NN.png (one per trial).

    if nargin < 2 || isempty(outSub), outSub = 'singleTrials'; end
    if nargin < 3 || isempty(epWin),  epWin = [-1.75 5.75]; end

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo, 'threshAnalysis'));
    L = labPaths();

    T = thresh_session_table(false); T = T(strcmp(T.sessID, sessID), :);
    if isempty(T), error('plot_session_singletrials_thresh:noSess', 'session not on disk: %s', sessID); end
    Sv = load(char(T.path(1))); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;

    labs   = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
    macIdx = find(cellfun(@(x) contains(x,'macBP'), labs));
    macLabs = labs(macIdx);
    if isempty(macIdx), error('plot_session_singletrials_thresh:noMacBP', '%s has no macBP channels', sessID); end
    bestMac = ''; if isfield(od,'bestMac'), bestMac = char(string(od.bestMac)); end
    biLocal = find(strcmp(macLabs, bestMac), 1);   % bestMac index within macLabs

    fs = od.fs; ts = round(od.TTL.start); nTr = numel(ts);
    s0 = round(epWin(1)*fs); s1 = round(epWin(2)*fs);

    % relative noise flags on bestMac (same rule as the pipeline)
    bestSig = double(od.data(macIdx(biLocal), :));
    NT = thresh_noise_trials(bestSig, fs, ts, epWin);     % uses calibrated K default
    times = NT.times;

    figDir = fullfile(L.figPath, sessID, 'threshTask', outSub);
    if ~isfolder(figDir), mkdir(figDir); end

    nSaved = 0;
    for t = 1:nTr
        if ~NT.ok(t), continue; end
        a = ts(t)+s0; b = ts(t)+s1;
        seg = double(od.data(macIdx, a:b));            % [nMac x nF]
        rng = max(seg,[],2) - min(seg,[],2);
        off = 1.15 * max([rng(:); eps]);

        f = figure('visible','off','position',[0 0 1150 720]); hold on;
        for c = 1:numel(macIdx)
            y = seg(c,:) + (c-1)*off;
            isBest = (c == biLocal);
            plot(times, y, 'color', [0.12 0.12 0.12], 'linewidth', 0.5 + 0.4*isBest);
            if isBest
                yr = y; yr(~NT.flagMask(:,t)') = NaN; plot(times, yr, 'r', 'linewidth', 1.1);
            end
            lbl = macLabs{c}; if isBest, lbl = [lbl ' (best)']; end %#ok<AGROW>
            if isBest, fw = 'bold'; lcol = [0.7 0 0]; else, fw = 'normal'; lcol = [0 0 0]; end
            text(times(1), (c-1)*off, [' ' lbl], 'fontsize', 9, ...
                'verticalalignment','bottom', 'interpreter','none', 'color', lcol, 'fontweight', fw);
        end
        xline(0, 'b', 'linewidth', 1);
        set(gca, 'ytick', []);
        ylim([min(seg(1,:))-0.3*off, (numel(macIdx)-1)*off + max(seg(end,:)) + 0.3*off]);
        xlim([times(1) times(end)]);
        xlabel('Time from TTL.start (ms)'); ylabel('bipolar channels (raw uV, stacked)');
        verdict = 'kept'; if NT.noisy(t), verdict = 'REJECTED'; end
        title(sprintf('%s   trial %d/%d   [%s]   (red = relative sharp deflection zd>%.1f on bestMac %s)', ...
            sessID, t, nTr, verdict, NT.K, bestMac), 'interpreter', 'none');
        saveas(f, fullfile(figDir, sprintf('trial_%02d.png', t))); close(f);
        nSaved = nSaved + 1;
    end
    fprintf('saved %d/%d per-trial plots to %s (rejected %d trials)\n', nSaved, nTr, figDir, sum(NT.noisy(NT.ok)));
end
