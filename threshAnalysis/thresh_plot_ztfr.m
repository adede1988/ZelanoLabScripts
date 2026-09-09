function thresh_plot_ztfr(map, timesMs, freqs, clim, resp, respT, dispWin, ttl, outPng)
% THRESH_PLOT_ZTFR  Spectrogram (log-freq, raw Hz labels, z heat) with the mean
%   respiration trace overlaid on a right-hand axis. Called once per odor
%   condition (low / med / air), all finalOnset-locked.

    m = timesMs >= dispWin(1) & timesMs <= dispWin(2);
    f = figure('visible','off','position',[0 0 880 540]);
    ax1 = axes('Parent', f);
    imagesc(ax1, timesMs(m), 1:numel(freqs), map(:,m)); axis(ax1, 'xy');

    tickHz = [2 4 8 16 32 64 120]; tickHz = tickHz(tickHz >= freqs(1) & tickHz <= freqs(end));
    yt = interp1(freqs, 1:numel(freqs), tickHz, 'linear', 'extrap');
    set(ax1, 'ytick', yt, 'yticklabel', string(tickHz));
    if ~isempty(clim) && all(isfinite(clim)) && clim(1) < clim(2), caxis(ax1, clim); end
    colormap(ax1, jet); cb = colorbar(ax1); ylabel(cb, 'Power (z vs baseline)');
    xlabel(ax1, 'Time (ms)'); ylabel(ax1, 'Frequency (Hz)');
    title(ax1, ttl, 'interpreter', 'none');
    xlim(ax1, dispWin); hold(ax1, 'on'); xline(ax1, 0, 'k', 'linewidth', 1.5);
    drawnow;

    % respiration overlay on a transparent axis sharing ax1's position
    if ~isempty(resp)
        mr = respT >= dispWin(1) & respT <= dispWin(2);
        pos = get(ax1, 'Position');
        ax2 = axes('Parent', f, 'Position', pos, 'Color', 'none', ...
                   'YAxisLocation', 'right', 'XTick', [], 'XLim', dispWin, ...
                   'XColor', 'none', 'YColor', [0 0 0]);
        hold(ax2, 'on');
        plot(ax2, respT(mr), resp(mr), 'color', [0 0 0], 'linewidth', 1.75);
        ylabel(ax2, 'mean respiration (a.u.)');
        ax2.XLim = dispWin;
    end
    saveas(f, outPng); close(f);
end
