function cue_plot_fooof(R, sessID, outPng)
% CUE_PLOT_FOOOF  Overlay the periodic (flattened) spectra of all macBP channels.
%   Periodic spectrum = 10*log10(PSD / aperiodic fit) = power above the 1/f.
%   bestMac is drawn bold; the 30-58 Hz gamma band is marked.

    f = figure('visible','off','position',[0 0 820 520]);
    hold on;
    nMac = numel(R.labels);
    cols = lines(max(nMac,1));
    for m = 1:nMac
        y = 10*log10(R.flattened(m,:));
        plot(R.freq, y, 'color', cols(m,:), 'linewidth', 1.2, ...
             'DisplayName', R.labels{m});
    end
    yb = 10*log10(R.flattened(R.bestIdx,:));
    plot(R.freq, yb, 'k', 'linewidth', 2.6, ...
         'DisplayName', sprintf('%s (best)', R.bestMac));

    yl = ylim;
    patch([R.gammaBand(1) R.gammaBand(2) R.gammaBand(2) R.gammaBand(1)], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
          'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
    set(gca,'children',flipud(get(gca,'children')));   % band behind lines

    xlim([2 120]); ylim(yl);
    xlabel('Frequency (Hz)'); ylabel('Periodic power (dB over aperiodic)');
    legend('show','location','northeast','interpreter','none');
    title(sprintf('%s  macBP periodic spectra  (best=%s, %s)', ...
        sessID, R.bestMac, R.selectionMethod), 'interpreter','none');
    saveas(f, outPng); close(f);
end
