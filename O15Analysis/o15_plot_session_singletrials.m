function o15_plot_session_singletrials(sessID, outSub, padPreS, padPostS)
% O15_PLOT_SESSION_SINGLETRIALS  Noise-diagnostic: one figure PER O15 TRIAL (the
%   ~15 odor trials) of a session's bipolar (macBP) channels, stacked, spanning
%   the whole trial (first sniff -> last/confirm sniff +/- pad). bestMac is drawn
%   bolder and the RELATIVE sharp-deflection windows are overlaid in RED on it,
%   so you can see exactly WHERE the detector fires against the raw signal. Each
%   sniff onset in the trial is marked by a vertical line colored by type
%   (start = green, free = blue, confirm = magenta); sniffs the rule REJECTS get
%   a red x. Ported from cueAnalysis/plot_session_singletrials (per trialStart
%   trial) to the O15 sniff structure.
%
%   The overlay uses the EXACT same relative rule as the pipeline: it calls
%   cue_noise_trials (with K = o15_noise_K) for the per-sniff reject flags, and
%   reconstructs the continuous per-sample flag from that call's own med/sigma/K
%   so the picture cannot drift from what run_o15_ztfr / run_o15_gamma_epochs do.
%
%   o15_plot_session_singletrials(sessID [, outSub, padPreS, padPostS])
%   Saves <figs>\<id>\O15\<outSub>\trial_NN.png (one per trial). Standalone
%   diagnostic; reads the final only, writes only these PNGs.

    if nargin < 2 || isempty(outSub),   outSub   = 'singleTrials'; end
    if nargin < 3 || isempty(padPreS),  padPreS  = 1.75; end
    if nargin < 4 || isempty(padPostS), padPostS = 3.0;  end

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'cueAnalysis')); addpath(fullfile(repo,'O15Analysis'));
    L = labPaths();

    T = o15_session_table(false); T = T(strcmp(T.sessID, sessID) & T.onDisk, :);
    if isempty(T), error('o15_plot_session_singletrials:noSess','session not on disk: %s', sessID); end
    Sv = load(char(T.path(1))); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;

    labs   = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
    macIdx = find(cellfun(@(x) contains(x,'macBP'), labs));
    macLabs = labs(macIdx);
    if isempty(macIdx), error('o15_plot_session_singletrials:noMacBP','%s has no macBP channels', sessID); end
    bestMac = ''; if isfield(od,'bestMac') && ~isempty(od.bestMac), bestMac = char(string(od.bestMac)); end
    biLocal = find(strcmp(macLabs, bestMac), 1);
    if isempty(biLocal)
        % Do NOT silently fall back to macBP1 -- that would diagnose the WRONG
        % channel (a flat channel inflates the robust-z and over-flags). A missing
        % bestMac means the final was re-preprocessed after FOOOF; re-run FOOOF.
        error('o15_plot_session_singletrials:noBestMac', ...
            ['%s has no usable bestMac (field missing or unmatched) -- the final was ' ...
             'likely re-preprocessed after FOOOF. Re-run run_o15_fooof_one(''%s'') first.'], ...
            sessID, sessID);
    end

    fs = od.fs; Tn = size(od.data,2);
    bd = od.behDat; sl = string(bd.sniffLabel); fo = round(bd.finalOnset); nn = bd.n;

    % per-sniff reject flags via the SAME pipeline rule (K from o15_noise_K)
    bestSig = double(od.data(macIdx(biLocal), :));
    NT = cue_noise_trials(bestSig, fs, fo, [], o15_noise_K());
    % continuous per-sample flag, reconstructed from NT's own med/sigma/win/K
    d  = movmax(bestSig, NT.win) - movmin(bestSig, NT.win);
    zd = (d - NT.med) ./ NT.sigma;
    flagFull = zd > NT.K;

    figDir = fullfile(L.figPath, sessID, 'O15', outSub);
    if ~isfolder(figDir), mkdir(figDir); end

    typeCol = containers.Map({'start','free','confirm'}, {[0 0.6 0],[0 0.2 0.9],[0.8 0 0.8]});
    utr = unique(nn(:))';
    nSaved = 0; totRej = 0;
    for u = utr
        rws = find(nn == u);
        if isempty(rws), continue; end
        ons = fo(rws);
        a = max(1,  round(min(ons) - padPreS*fs));
        b = min(Tn, round(max(ons) + padPostS*fs));
        if b - a < round(0.5*fs), continue; end
        seg = double(od.data(macIdx, a:b));            % [nMac x nW]
        times = ((a:b) - min(ons)) / fs * 1000;        % ms rel. to first sniff in trial
        rng = max(seg,[],2) - min(seg,[],2);
        off = 1.15 * max([rng(:); eps]);

        f = figure('visible','off','position',[0 0 1300 760]); hold on;
        for c = 1:numel(macIdx)
            y = seg(c,:) + (c-1)*off;
            isBest = (c == biLocal);
            plot(times, y, 'color', [0.12 0.12 0.12], 'linewidth', 0.5 + 0.4*isBest);
            if isBest
                yr = y; fm = flagFull(a:b); yr(~fm) = NaN;
                plot(times, yr, 'r', 'linewidth', 1.2);    % relative sharp-deflection windows
            end
            lbl = macLabs{c}; fw = 'normal'; lc = [0 0 0];
            if isBest, lbl = [lbl ' (best)']; fw = 'bold'; lc = [0.7 0 0]; end %#ok<AGROW>
            text(times(1), (c-1)*off, [' ' lbl], 'fontsize',9, 'verticalalignment','bottom', ...
                'interpreter','none', 'color',lc, 'fontweight',fw);
        end

        yl = [min(seg(1,:)) - 0.3*off, (numel(macIdx)-1)*off + max(seg(end,:)) + 0.6*off];
        nRej = 0;
        for r = rws'
            xo = (fo(r) - min(ons)) / fs * 1000;
            tp = char(sl(r)); col = [0.5 0.5 0.5]; if isKey(typeCol,tp), col = typeCol(tp); end
            plot([xo xo], yl, '-', 'color', col, 'linewidth', 1.0);
            if NT.noisy(r)
                nRej = nRej + 1;
                plot(xo, yl(2), 'rx', 'markersize', 11, 'linewidth', 1.8);
            end
        end
        totRej = totRej + nRej;
        set(gca,'ytick',[]); ylim(yl); xlim([times(1) times(end)]);
        xlabel('Time from first sniff in trial (ms)'); ylabel('macBP channels (raw uV, stacked)');
        nS = numel(rws);
        nSt = sum(sl(rws)=="start"); nFr = sum(sl(rws)=="free"); nCo = sum(sl(rws)=="confirm");
        title(sprintf(['%s   trial %d   sniffs %d (start %d/free %d/confirm %d), REJECTED %d   ' ...
            '[red = relative zd>%.1f on bestMac %s; lines: green start / blue free / magenta confirm; x = rejected]'], ...
            sessID, u, nS, nSt, nFr, nCo, nRej, NT.K, bestMac), 'interpreter','none','fontsize',8);
        saveas(f, fullfile(figDir, sprintf('trial_%02d.png', u))); close(f);
        nSaved = nSaved + 1;
    end
    fprintf('%s: saved %d/%d per-trial diagnostics to %s (%d sniffs rejected total, K=%.1f)\n', ...
        sessID, nSaved, numel(utr), figDir, totRej, NT.K);
end
