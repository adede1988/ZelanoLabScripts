function chirp_fig_trial(tfr, R, ch, i, sessID, chi, figDir, C, chOut)
% CHIRP_FIG_TRIAL  Single-trial FASLT spectrogram + ridge + anchors + transition + chirplet atoms.
%   chirp_fig_trial(tfr, R, ch, i, sessID, chi, figDir, C, chOut)
%   Saves sub-<sessID>_ch-<chi>_trial-<NNN>.png under figDir (spec 2.5/6.3). chOut optional:
%   per-trial chirplet output (to overlay the dominant atom trajectory).

    if ~isfolder(figDir), mkdir(figDir); end
    f = tfr.freqs; t = tfr.tMs; P = tfr.power(:,:,i);
    if all(~isfinite(P(:))), return; end

    fh = figure('Visible','off','Position',[100 100 760 460]);
    imagesc(t, f, 10*log10(P + eps)); axis xy; hold on;
    colormap(parula); cb = colorbar; cb.Label.String = 'power (dB)';
    ylim(C.faslt.range); xlim([t(1) t(end)]);
    xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');

    T = R.trial(i);
    if ~isempty(T.fhat), plot(t, T.fhat, 'w-', 'LineWidth', 1.6); end
    if isfield(T,'fhat2') && ~isempty(T.fhat2) && T.nRidgesSig>1
        plot(t, T.fhat2, 'w--', 'LineWidth', 1.0);
    end
    yline(R.anchors.f_hi, 'r:', 'f_{hi}', 'LineWidth', 1);
    yline(R.anchors.f_lo, 'r:', 'f_{lo}', 'LineWidth', 1);
    if ~isempty(T.transitionWin)
        xline(T.transitionWin(1), 'g-', 'LineWidth', 1);
        xline(T.transitionWin(2), 'g-', 'LineWidth', 1);
    end
    if nargin>=9 && ~isempty(chOut) && isstruct(chOut) && ~isempty(chOut.atoms)
        A = chOut.atoms;   % [tc_s fc_Hz chirp sigma amp phase eFrac]
        plot(A(:,1)*1000, A(:,2), 'mo-', 'MarkerFaceColor','m', 'MarkerSize',4, 'LineWidth',1);
    end
    cls = ''; if nargin>=9 && ~isempty(chOut), cls = chOut.classification; end
    title(sprintf('%s  ch=%s  trial %d   nRidgesSig=%d  dip=%.2f  %s', ...
        strrep(sessID,'_','\_'), strrep(ch,'_','\_'), i, T.nRidgesSig, T.powerDipStat, cls), ...
        'Interpreter','tex');

    fn = fullfile(figDir, sprintf('sub-%s_ch-%d_trial-%03d.png', sessID, chi, i));
    exportgraphics(fh, fn, 'Resolution', 110);
    close(fh);
end
